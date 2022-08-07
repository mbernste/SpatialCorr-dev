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
import time

import scHOT

def main():
    usage = "%prog [options] input_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="Output file")
    (options, args) = parser.parse_args()

    gene_1, gene_2 = args[0].split(',')

    out_file = options.out_file

    print("Loading data...")
    print("Loading Dino-normalized data...")
    df_dino = pd.read_csv('../SpatialLIBD_data/spatiallbd_dino.tsv.gz', sep='\t', index_col=0)
    df_dino = df_dino.transpose()
    df_dino.index = [x.replace('.', '-') for x in df_dino.index]
    adata = AnnData(df_dino)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    print(adata.X.shape)

    # Load metadata
    meta_df = pd.read_csv('../SpatialLIBD_data/spatialLIBD_coldata.tsv', sep='\t')
    meta_df = meta_df.set_index('barcode')
    meta_df = meta_df.loc[meta_df['sample_name'] == 151507]
    meta_df = meta_df.loc[adata.obs.index]
    meta_df = meta_df.loc[meta_df['tissue'] == 1]
    meta_df['layer_guess'] = [str(x) for x in meta_df['layer_guess']]
    assert tuple(meta_df.index) == tuple(adata.obs.index)

    # Run LLRT on each pair of most variable genes
    start = time.time()
    p_val = run_scHOT(gene_1, gene_2, meta_df, adata, adata.X)
    end = time.time()
    elapsed_time = end - start
    print("Elapsed time: ", elapsed_time)
    with open(out_file, 'w') as f:
        json.dump(
            {
                'elapsed_time': elapsed_time,
                'p_val': p_val
            },
            f,
            indent=True
        )

def run_scHOT(gene_1, gene_2, meta_df, ad, X):
    ind_1 = list(ad.var.index).index(gene_1)
    ind_2 = list(ad.var.index).index(gene_2)
    expr_1 = X[:,ind_1]
    expr_2 = X[:,ind_2]
    df_schot = pd.DataFrame(
        data={
            gene_1: expr_1,
            gene_2: expr_2,
            'imagerow': meta_df['imagerow'],
            'imagecol': meta_df['imagecol']
        },
        index=meta_df.index
    )
    df_schot_res = scHOT.run(
        df_schot, 
        [gene_1, gene_2], 
        f'{gene_1}_{gene_2}', 
        tmp_dir='./scHOT_tmp', 
        script_dir='.', 
        verbose=True
    )
    p_val = df_schot_res.loc[gene_1]['pvalEstimated']
    print('scHOT p-value: ', p_val)
    return p_val 

if __name__ == '__main__':
    main()
