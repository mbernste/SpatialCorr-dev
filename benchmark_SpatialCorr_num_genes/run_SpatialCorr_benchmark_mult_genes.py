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

#sys.path.append('..')
sys.path.append('../../spatialcorr')

import spatialcorr

N_PROCS = 25
GENE_SET_SIZE_MAX = 10
GENE_SET_SIZE_MIN = 5

def main():
    usage = "%prog [options] input_file"
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_file", help="Output file")
    parser.add_option("-n", "--no_mc_pvals", action='store_true', help="Don't use Monte Carlo p-values")
    parser.add_option("-b", "--run_br_test", action='store_true', help="Run the BR-test")
    (options, args) = parser.parse_args()

    mc_pvals = not options.no_mc_pvals
    br_test = options.run_br_test
    out_file = options.out_file

    print("Loading top genes...")
    top_genes = list(pd.read_csv('top_genes.tsv', sep='\t')['gene'])

    print("Loading data...")
    print("Loading Dino-normalized data...")
    df_dino = pd.read_csv('./SpatialLIBD_data/spatiallbd_dino.tsv.gz', sep='\t', index_col=0)
    df_dino = df_dino.transpose()
    df_dino.index = [x.replace('.', '-') for x in df_dino.index]
    adata = AnnData(df_dino)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)
    print(adata.X.shape)

    # Load metadata
    meta_df = pd.read_csv('./SpatialLIBD_data/spatialLIBD_coldata.tsv', sep='\t')
    meta_df = meta_df.set_index('barcode')
    meta_df = meta_df.loc[meta_df['sample_name'] == 151507]
    meta_df = meta_df.loc[adata.obs.index]
    meta_df = meta_df.loc[meta_df['tissue'] == 1]
    meta_df['layer_guess'] = [str(x) for x in meta_df['layer_guess']]
    assert tuple(meta_df.index) == tuple(adata.obs.index)

    # Add metadata to AnnData dataset
    adata.obs['row'] = meta_df['row']
    adata.obs['col'] = meta_df['col']
    adata.obs['imagerow'] = meta_df['imagerow']
    adata.obs['imagecol'] = meta_df['imagecol']
    adata.obs['layer_guess'] = meta_df['layer_guess']

    # The gene to use (it doesn't effect the experiment)
    GENE = 'SNAP25'
    NUM_GENES = [2, 4, 8, 16, 32, 64, 128]

    all_elapsed = []
    for num_genes in NUM_GENES:
        print(f"Running with {num_genes} genes...")
        #genes = [GENE for i in range(num_genes)]

        genes = top_genes[:num_genes]

        print(genes)

        start = time.time()
        p_val, additional = spatialcorr.run_test(
            adata,
            genes,
            5,
            cond_key='layer_guess',
            contrib_thresh=10,
            row_key='row',
            col_key='col',
            verbose=1,
            n_procs=1,
            max_perms=10,
            compute_spotwise_pvals=False,
            mc_pvals=False,
            run_bhr=False
        )
        end = time.time()
 
        elapsed_time = end - start
        print("Elapsed: ", elapsed_time)
        all_elapsed.append(elapsed_time)


    df_result = pd.DataFrame(
        data={
            'number_of_genes': NUM_GENES,
            'elapsed_time': all_elapsed
        }
    ).set_index('number_of_genes')
    df_result.to_csv('results.tsv', sep='\t')

if __name__ == '__main__':
    main()
