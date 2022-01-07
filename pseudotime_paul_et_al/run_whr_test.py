import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy.io import mmread
import scanpy as sc
from anndata import AnnData
import importlib
import sys
from collections import defaultdict
from optparse import OptionParser
import json
import gc

#sys.path.append('..')
sys.path.append('../../spatialcorr')

import spatialcorr
#import normalize_deviance_residuals as ndr

N_PROCS = 25
GENE_SET_SIZE_MAX = 15
GENE_SET_SIZE_MIN = 5

def main():
    usage = "%prog [options] input_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--between_cell_types", action="store_true", help="Test for differences in correlation between cell types.")
    parser.add_option("-g", "--gene_sets_group_file", help="File storing set of gene sets to run experiment on")
    parser.add_option("-i", "--high_expression", action="store_true", help="Test only sets of genes with high expression")
    parser.add_option("-p", "--num_procs", help="Number of processes")
    parser.add_option("-j", "--override_sig_joint", action='store_true', help="Run pairwise tests regardless of joint p-value")
    parser.add_option("-s", "--sigma", help="Sigma")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    between_cell_types = options.between_cell_types
    only_high = options.high_expression
    override_sig_joint = options.override_sig_joint
    gene_sets_f = options.gene_sets_group_file
    out_file = options.out_file
    sigma = float(options.sigma)

    print(f"Using sigma={sigma}")

    if options.num_procs is not None:
        n_procs = int(options.num_procs)
    else:
        n_procs = N_PROCS


    # Load cell-cell distances
    print("Loading distance matrix...")
    df_dists = pd.read_csv('./Paul15_dist_matrix.tsv', sep='\t', index_col=0)
    dist_matrix = np.array(df_dists)

    # Load Dino data
    print("Loading data Dino-normalized data...")
    df_dino = pd.read_csv('./Paul15_dino.tsv', sep='\t', index_col=0)
    df_dino = df_dino.transpose()
    adata = AnnData(df_dino)
    assert tuple(df_dino.index) == tuple(df_dists.index)

    # Load counts
    print("Loading data counts data...")
    df_counts = pd.read_csv('./Paul15_counts.tsv', sep='\t', index_col=0)
    df_counts = df_counts.transpose()
    adata_counts = AnnData(df_counts)
    assert tuple(df_counts.index) == tuple(df_dists.index)

    # Load metadata
    print("Loading metadata...")
    df_meta = pd.read_csv('./Paul15_meta.tsv', sep='\t', index_col=0)
    adata.obs['cluster'] = df_meta['cluster']
    assert tuple(df_meta.index) == tuple(df_dists.index)

    # Get all genes expression in > 0.2 of spots
    print(f"Computing percentage of spots expressing each gene...")
    gene_totals = np.sum((np.array(adata_counts.X) > 0).astype(int), axis=0)
    gene_totals = gene_totals / adata_counts.X.shape[0]
    genes_gt_0_2_spots = [
        (gene, total)
        for gene, total in zip(adata_counts.var.index, gene_totals)
        if total > 0.2
    ]
    gene_to_total = {
        gene: total
        for gene, total in zip(adata_counts.var.index, gene_totals)
    }
    print(f"{len(genes_gt_0_2_spots)} candidate genes expr > 0.2 spots...")

    # Load gene sets
    gene_set_to_genes_GO_bp = _parse_gene_sets('./gene_sets/c5.bp.v7.1.symbols.gmt')

    # Garbage collect
    adata_counts = None
    gc.collect()

    # Load gene sets
    gene_sets = list(pd.read_csv(gene_sets_f, header=None, sep='\t', index_col=0).index)
    gene_sets = [x.strip() for x in gene_sets]
    print(f"{len(gene_sets)} gene sets: {gene_sets}")

    # Compute kernel matrix
    kernel_matrix = spatialcorr.statistical_test._compute_kernel_matrix(
        adata.obs,
        sigma,
        condition_on_cell_type=True,
        cell_type_key='cluster',
        dist_matrix=dist_matrix
    )

    # Map lower-case version of genes to original
    lower_to_orig = {
        g.lower(): g
        for g in adata.var.index
    }

    gene_set_to_genes = {}
    gene_set_to_pval = {}
    gene_set_to_ct_to_pval = {}
    gene_set_to_pairwise_pvals = defaultdict(lambda: {})
    gene_set_to_pairwise_ct_to_pvals = defaultdict(lambda: {})
    cache_genes_to_pval = {}
    for gs_i, gs in enumerate(gene_sets):
        top_expr = sorted(
            [
                lower_to_orig[gene.lower()]
                for gene in gene_set_to_genes_GO_bp[gs]
                if gene.lower() in lower_to_orig
            ], 
            key=lambda x: gene_to_total[x], 
            reverse=True
        )
        # Grab most expressed genes
        if len(top_expr) > GENE_SET_SIZE_MAX:
            top_expr = top_expr[:GENE_SET_SIZE_MAX]
        # Filter genes expressed in at least 0.25 spots
        genes_meet_thresh = sorted([
            gene
            for gene in top_expr
            if gene_to_total[gene] > 0.2
        ])
        # Skip if we have too few genes
        if len(genes_meet_thresh) < GENE_SET_SIZE_MIN:
            print(f"Gene set {gs}. Only {len(genes_meet_thresh)} genes in the gene set. Skipping gene set...")
            continue 

        print(f"Gene set {gs}. Genes meeting threshold: {genes_meet_thresh}")
        gene_set_to_genes[gs] = genes_meet_thresh
        inds = [
            list(adata.var.index).index(gene)
            for gene in genes_meet_thresh
        ]
        p_val, additional = spatialcorr.run_test(
            adata,
            genes_meet_thresh,
            sigma,
            cond_key='cluster',
            precomputed_kernel=kernel_matrix,
            contrib_thresh=10,
            verbose=1,
            n_procs=n_procs,
            max_perms=10000, # TODO CHANGE TO 10000
            compute_spotwise_pvals=True,
            run_bhr=False
        )
        ct_to_p_val = additional['region_to_p_val']
        print("p-value: ", p_val)
     
        gene_set_to_pval[gs] = p_val
        gene_set_to_ct_to_pval[gs] = ct_to_p_val        

        # Run pairwise test
        if p_val < 0.05 or override_sig_joint:
            for g1_i, g1 in enumerate(sorted(genes_meet_thresh)):
                for g2_i, g2 in enumerate(sorted(genes_meet_thresh)):
                    if g1_i >= g2_i:
                        continue   
                    print("Running test on ", g1, g2)
                    inds = [
                        list(adata.var.index).index(gene)
                        for gene in [g1, g2]
                    ]

                    # Check the cache to see if we've already computed this p-value
                    if frozenset({g1, g2}) in cache_genes_to_pval.keys():
                        pair_p_val = cache_genes_to_pval[frozenset({g1, g2})]
                        print(f"Already computed a p-value for {g1} and {g2}: {pair_p_val}")
                    else:
                        pair_p_val, pair_additional = spatialcorr.run_test(
                            adata,
                            [g1, g2],
                            sigma,
                            cond_key='cluster',
                            precomputed_kernel=kernel_matrix,
                            contrib_thresh=10,
                            verbose=1,
                            n_procs=n_procs,
                            max_perms=100, # TODO CHANGE TO 1000
                            compute_spotwise_pvals=True,
                            run_bhr=False
                        )
                        pair_ct_to_p_val = pair_additional['region_to_p_val']
                        cache_genes_to_pval[frozenset({g1, g2})] = pair_p_val # Record the new p-value
                 
                    print("p-value: ", pair_p_val)
                    gene_set_to_pairwise_pvals[gs][f'{g1},{g2}'] = pair_p_val
                    gene_set_to_pairwise_ct_to_pvals[gs][f'{g1},{g2}'] = pair_ct_to_p_val

    # Write output
    with open(out_file, 'w') as f:
        json.dump(
            {
                'gene_set_to_genes': gene_set_to_genes,
                'gene_set_to_pval': gene_set_to_pval,
                'gene_set_to_pairwise_pvals': gene_set_to_pairwise_pvals,
                'gene_set_to_pairwise_clust_to_pval': gene_set_to_pairwise_ct_to_pvals
            },
            f,
            indent=True
        )


def _parse_gene_sets(gene_sets_f):
    gene_set_to_genes = {}
    with open(gene_sets_f, 'r') as f:
        for l in f:
            toks = l.split('\t')
            gene_set = toks[0]
            genes = [x.strip() for x in toks[2:]]
            gene_set_to_genes[gene_set] = genes
    return gene_set_to_genes

if __name__ == '__main__':
    main()
