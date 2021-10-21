import os
from os.path import join
import pandas as pd
from optparse import OptionParser
import sys
import numpy as np
import scanpy as sc
from anndata import AnnData

import spatialcorr

def main():
    parser = OptionParser()
    parser.add_option("-o", "--out_file", help="Output file")
    parser.add_option("-c", "--condition_column", help="Column in input datasets on which to condition")
    parser.add_option("-w", "--between_conditions", action='store_true', help="If true, run between cell type test")
    parser.add_option("-b", "--bandwidth", help="The bandwidth of the kernel to use")
    parser.add_option("-r", "--regress_size_factors", action='store_true', help="Regress out size-factors from deviance residuals")
    parser.add_option("-t", "--contrib_thresh", help="Spotwise effective number of data points threshold")
    parser.add_option("-v", "--std_var", action='store_true', help="Standardize variance")
    parser.add_option("-f", "--feats", help="Features to use")
    parser.add_option("-n", "--n_procs", help="Number of processes")
    (options, args) = parser.parse_args()

    in_dir = args[0]
    sigma = float(options.bandwidth)
    n_procs = int(options.n_procs)
    condition_col = options.condition_column
    between_conds = options.between_conditions
    std_var = options.std_var
    regress_size_factors = options.regress_size_factors
    feats = options.feats
    out_f = options.out_file

    if not options.contrib_thresh:
        contrib_thresh = 30
    else:
        contrib_thresh = float(options.contrib_thresh)

    # Whether to condition on cell type or some other covariate
    condition = not (condition_col is None)

    print()
    print("Running test with the following parameters:")
    print(f"Bandwidth: {sigma}")
    print(f"Testing for correlation between conditions: {between_conds}")
    print(f"Features used: {feats}")
    print(f"Condition on covariate: {condition}")
    print(f"Covariate conditioned on: {condition_col}")
    print(f"Number of processes: {n_procs}")
    print(f"Regress out size-factors: {regress_size_factors}")
    print(f"Spotwise effective contribution threshold: {contrib_thresh}")
    print(f"Standardize the variance: {std_var}")
    print()

    if feats == 'dino':
        feat_1 = 'dino_gene_1'
        feat_2 = 'dino_gene_2'
    elif feats == 'dev_res':
        feat_1 = 'dev_res_gene_1'
        feat_2 = 'dev_res_gene_2'

    da = []
    for elem_i, elem in enumerate(os.listdir(in_dir)):
        in_file = join(in_dir, elem)
        print(f"Running test on data in {in_file}...")

        # Read simulated data file
        df_sim = pd.read_csv(in_file, sep='\t', index_col=0)

        expr = np.array(df_sim[[
            feat_1, 
            feat_2
        ]])
        adata = AnnData(
            X=expr,
            obs=df_sim[['row', 'col']],
            var=pd.DataFrame(
                data=['gene_1', 'gene_2'],
                index=['gene_1', 'gene_2']
            )
        )
        adata.obs['size_factor'] = df_sim['size_factors']

        if 'cluster' in set(df_sim.columns):
            adata.obs['cluster'] = df_sim['cluster']        

        if regress_size_factors:
            print("Regressing out size-factors...")
            sc.pp.regress_out(adata, 'size_factor')

        p_val, t_obs, t_nulls, keep_inds, spotwise_t_obs, spotwise_t_nulls, spot_p_vals = spatialcorr.run_test(
            adata,
            ['gene_1', 'gene_2'],
            sigma,
            cond_key=condition_col,
            contrib_thresh=contrib_thresh,
            row_key='row',
            col_key='col',
            verbose=1,
            n_procs=n_procs,
            test_between_conds=between_conds,
            standardize_var=std_var
        ) 
        print(f"p-value: {p_val}")

        elem_name = elem.split('.')[0]
        da.append([elem_name, p_val])

    df_out = pd.DataFrame(
        data=da,
        columns=['dataset', 'p_value']
    )
    df_out = df_out.set_index('dataset')
    df_out.to_csv(out_f, sep='\t')

if __name__ == '__main__':
    main()
