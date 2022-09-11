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

#sys.path.append('..')
sys.path.append('../../spatialcorr')

import spatialcorr

N_PROCS = 25
GENE_SET_SIZE_MAX = 10
GENE_SET_SIZE_MIN = 5

def main():
    usage = "%prog [options] input_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--condition_cluster", action="store_true", help="Condition on cluster.")
    parser.add_option("-b", "--between_cell_types", action="store_true", help="Test for differences in correlation between cell types.")
    parser.add_option("-g", "--gene_sets_group_file", help="File storing set of gene sets to run experiment on")
    parser.add_option("-j", "--override_sig_joint", action='store_true', help="Run pairwise tests regardless of joint p-value")
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    print("Loading Dino-normalized data...")
    df_dino = pd.read_csv('../SpatialLIBD_data/spatiallbd_dino.tsv.gz', sep='\t', index_col=0)
    df_dino = df_dino.transpose()
    df_dino.index = [x.replace('.', '-') for x in df_dino.index]
    adata = AnnData(df_dino)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    print(adata.X.shape)

    print("Loading counts data...")
    adata_raw = sc.read_10x_h5('../SpatialLIBD_data/151507_filtered_feature_bc_matrix.h5')
    adata_raw.X = np.array(adata_raw.X.todense())

    # Load metadata
    meta_df = pd.read_csv('../SpatialLIBD_data/spatialLIBD_coldata.tsv', sep='\t')
    meta_df = meta_df.set_index('barcode')
    meta_df = meta_df.loc[meta_df['sample_name'] == 151507]
    meta_df = meta_df.loc[adata.obs.index]
    meta_df = meta_df.loc[meta_df['tissue'] == 1]
    meta_df['layer_guess'] = [str(x) for x in meta_df['layer_guess']]
    assert tuple(meta_df.index) == tuple(adata.obs.index)

    adata.obs['row'] = meta_df['row']
    adata.obs['col'] = meta_df['col']
    adata.obs['imagerow'] = meta_df['imagerow']
    adata.obs['imagecol'] = meta_df['imagecol']
    adata.obs['layer_guess'] = meta_df['layer_guess']

    # Get all genes expression in > 0.5 of spots
    print("Computing all genes expressed in at least half of all spots...")
    #gene_expr_total = np.sum(adata.X, axis=0)
    gene_totals = np.sum((adata_raw.X > 0).astype(int), axis=0)
    gene_totals = gene_totals / adata_raw.X.shape[0]
    genes_gt_0_2_spots = [
        (gene, total)
        for gene, total in zip(adata_raw.var.index, gene_totals)
        if total > 0.2
    ]

    gene_to_total = {
        gene: total
        for gene, total in zip(adata_raw.var.index, gene_totals)
    }
 
    dino_genes = set(adata.var.index)
    top_expr_genes = sorted(
        [
            x[0]
            for x in genes_gt_0_2_spots
            if x[0] in dino_genes
        ],
        key=lambda x: gene_to_total[x],
        reverse=True
    )

    da = [
        (gene, gene_to_total[gene])
        for gene in top_expr_genes
    ]
    df = pd.DataFrame(
        data=da,
        columns=['gene', 'fraction_spots_detected']
    )
    df = df.set_index('gene')
    df.to_csv('top_genes.tsv', sep='\t')

if __name__ == '__main__':
    main()
