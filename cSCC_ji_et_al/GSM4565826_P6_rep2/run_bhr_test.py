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
    parser.add_option("-c", "--condition_cluster", action="store_true", help="Condition on cluster.")
    parser.add_option("-b", "--between_cell_types", action="store_true", help="Test for differences in correlation between cell types.")
    parser.add_option("-g", "--gene_sets_group_file", help="File storing set of gene sets to run experiment on")
    parser.add_option("-n", "--num_clusts", help="Number of clusters to use")
    parser.add_option("-i", "--high_expression", action="store_true", help="Test only sets of genes with high expression")
    parser.add_option("-p", "--num_procs", help="Number of processes")
    parser.add_option("-j", "--override_sig_joint", action='store_true', help="Run pairwise tests regardless of joint p-value")
    parser.add_option("-s", "--sigma", action='store_true', help="Sigma (default: 5)")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    between_cell_types = options.between_cell_types
    condition_cluster = options.condition_cluster
    n_clust = options.num_clusts
    only_high = options.high_expression
    override_sig_joint = options.override_sig_joint
    gene_sets_f = options.gene_sets_group_file
    out_file = options.out_file

    if options.sigma is not None:
        sigma = float(options.sigma)
    else:
        sigma = 5

    if options.num_procs is not None:
        n_procs = int(options.num_procs)
    else:
        n_procs = N_PROCS

    # Load spatial metadata
    meta_df = pd.read_csv('../GSM4565826_P6_rep2_data/spatial/tissue_positions_list.csv')
    meta_df.columns = ['barcode', 'tissue', 'row', 'col', 'imgrow', 'imgcol']
    meta_df = meta_df.set_index('barcode')
    meta_df = meta_df.loc[meta_df['tissue'] == 1]

    # Load BayesSpace clusters
    df_bayes = pd.read_csv(f'../GSM4565826_P6_rep2_data/BayesSpace/BayesSpace_output.{n_clust}_clusts.tsv', sep='\t', index_col=0)
    df_bayes = df_bayes.loc[meta_df.index]
    meta_df[f'BayesSpace_{n_clust}'] = df_bayes['spatial.cluster']

    # Load Dino data
    print("Loading data Dino-normalized data...")
    df_dino = pd.read_csv('../GSM4565826_P6_rep2_data/dino/GSM4565826_P6_dino.tsv.gz', sep='\t', index_col=0)
    df_dino = df_dino.transpose()
    df_dino.index = [x.replace('.', '-') for x in df_dino.index]
    df_dino = df_dino.loc[meta_df.index]
    adata = AnnData(df_dino)
    assert tuple(df_dino.index) == tuple(adata.obs.index)
    assert tuple(df_dino.index) == tuple(meta_df.index)
    adata.obs = meta_df

    # Load counts data. Used for filtering
    print("Loading counts data...")
    adata_counts = sc.read_10x_mtx('../GSM4565826_P6_rep2_data/filtered_feature_bc_matrix')
    adata_counts.X = adata_counts.X.todense()
    assert adata_counts.shape[0] == adata.shape[0]

    # Get all genes expression in > 0.5 of spots
    print(f"Computing percentage of spots expressing each gene...")
    #print(np.array(adata_counts.X))
    gene_totals = np.sum((np.array(adata_counts.X) > 0).astype(int), axis=0)
    gene_totals = gene_totals / adata_counts.X.shape[0]

    # Get all genes expression in > 0.2 of spots
    print(f"Computing percentage of spots expressing each gene...")
    #print(np.array(adata_counts.X))
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

    gene_set_to_genes_GO_bp = _parse_gene_sets('./gene_sets/c5.bp.v7.1.symbols.gmt')

    # Garbage collect
    adata_counts = None
    gc.collect()

    # Load gene sets
    gene_sets = list(pd.read_csv(gene_sets_f, header=None, sep='\t', index_col=0).index)
    gene_sets = [x.strip() for x in gene_sets]
    print(f"{len(gene_sets)} gene sets: {gene_sets}")

    gene_set_to_genes = {}
    gene_set_to_pval = {}
    gene_set_to_pairwise_pvals = defaultdict(lambda: {})
    gene_set_to_keep_clusts = {}
    cache_genes_to_pval = {}
    for gs_i, gs in enumerate(gene_sets):
        top_expr = sorted(
            [
                gene 
                for gene in gene_set_to_genes_GO_bp[gs]
                if gene in adata.var.index
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
        if condition_cluster:
            cond_key=f'BayesSpace_{n_clust}'
        else:
            cond_key=None
        p_val, t_obs, t_nulls, kept_inds, obs_spot_lls, spotwise_t_nulls, ct_to_p_val = spatialcorr.run_test(
            adata,
            genes_meet_thresh,
            sigma,
            cond_key=cond_key,
            contrib_thresh=10,
            row_key='row',
            col_key='col',
            verbose=1,
            max_perms = 200,
            n_procs=n_procs,
            compute_spotwise_pvals=True,
            run_bhr=False
        )
        print("WHR p-value: ", p_val)

        keep_cts = set([
            ct
            for ct, pval in ct_to_p_val.items()
            if pval > 0.05
        ])
        print(f"Kept {len(keep_cts)} clusters after WHR test...")
        if len(keep_cts) < 2:
            continue       

        print("Filtering adata...")
        keep_ct_inds = [
            ind
            for ind, ct in zip(adata.obs.index, adata.obs[cond_key])
            if ct in keep_cts
        ]
        adata_bhr = adata[keep_ct_inds]

        p_val, t_obs, t_nulls, kept_inds, obs_spot_lls, spotwise_t_nulls, _ = spatialcorr.run_test(
            adata_bhr,
            genes_meet_thresh,
            sigma,
            cond_key=cond_key,
            contrib_thresh=10,
            row_key='row',
            col_key='col',
            verbose=1,
            n_procs=n_procs,
            #max_perms = 200, # TODO REMOVE
            run_bhr=True
        )
        print("BHR p-value: ", p_val)
 
        gene_set_to_pval[gs] = p_val
        gene_set_to_keep_clusts[gs] = sorted(keep_cts) 


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
                    if condition_cluster:
                        cond_key=f'BayesSpace_{n_clust}'
                    else:
                        cond_key=None

                    # Check the cache to see if we've already computed this p-value
                    if frozenset({g1, g2}) in cache_genes_to_pval.keys():
                        pair_p_val = cache_genes_to_pval[frozenset({g1, g2})]
                        print(f"Already computed a p-value for {g1} and {g2}: {pair_p_val}")
                    else:
                        pair_p_val, t_obs, t_nulls, kept_inds, obs_spot_lls, spotwise_t_nulls, _ = spatialcorr.run_test(
                            adata_bhr,
                            [g1, g2],
                            sigma,
                            cond_key=cond_key,
                            contrib_thresh=10,
                            row_key='row',
                            col_key='col',
                            verbose=1,
                            n_procs=n_procs,
                            max_perms=1000,
                            #max_perms=200, # TODO REMOVE
                            run_bhr=between_cell_types
                        )
                        cache_genes_to_pval[frozenset({g1, g2})] = pair_p_val # Record the new p-value
                 
                    print("p-value: ", pair_p_val)
                    gene_set_to_pairwise_pvals[gs][f'{g1},{g2}'] = pair_p_val

    # Write output
    with open(out_file, 'w') as f:
        json.dump(
            {
                'gene_set_to_genes': gene_set_to_genes,
                'gene_set_to_pval': gene_set_to_pval,
                'gene_set_to_pairwise_pvals': gene_set_to_pairwise_pvals,
                'gene_set_to_keep_clusts': gene_set_to_keep_clusts
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
