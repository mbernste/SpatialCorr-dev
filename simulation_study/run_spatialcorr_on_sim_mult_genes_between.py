import os
from os.path import join
import pandas as pd
from optparse import OptionParser
import sys
import numpy as np
import scanpy as sc
from anndata import AnnData
from collections import defaultdict
import json

import spatialcorr

def main():
    parser = OptionParser()
    parser.add_option("-o", "--out_file", help="Output file")
    parser.add_option("-c", "--condition_column", help="Column in input datasets on which to condition")
    parser.add_option("-w", "--between_conditions", action='store_true', help="If true, run between cell type test")
    parser.add_option("-b", "--bandwidth", help="The bandwidth of the kernel to use")
    parser.add_option("-t", "--contrib_thresh", help="Spotwise effective number of data points threshold")
    parser.add_option("-f", "--feats", help="Features to use")
    parser.add_option("-m", "--n_genes", help="Number of genes to run the test on")
    parser.add_option("-n", "--n_procs", help="Number of processes")
    (options, args) = parser.parse_args()

    in_dir = args[0]
    sigma = float(options.bandwidth)
    n_procs = int(options.n_procs)
    n_genes = int(options.n_genes)
    condition_col = options.condition_column
    between_conds = options.between_conditions
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
    print(f"Spotwise effective contribution threshold: {contrib_thresh}")
    print(f"Number of genes: {n_genes}")
    print()

    if feats == 'dino':
        all_feats = [
            f'dino_gene_{i}'
            for i in range(n_genes)
        ]
    elif feats == 'dev_res':
        all_feats = [
            f'dev_res_gene_{i}'
            for i in range(n_genes)
        ]
    gene_names = [f'gene_{i}' for i in range(n_genes)]

    da = []
    dataset_to_results = defaultdict(lambda: {})
    for elem_i, elem in enumerate(os.listdir(in_dir)):

        in_file = join(in_dir, elem)
        elem_name = elem.split('.')[0]

        print(f"Running test on data in {in_file}...")

        # Read simulated data file
        df_sim = pd.read_csv(in_file, sep='\t', index_col=0)

        expr = np.array(df_sim[all_feats])
        adata = AnnData(
            X=expr,
            obs=df_sim[['row', 'col']],
            var=pd.DataFrame(
                data=gene_names,
                index=gene_names
            )
        )
        adata.obs['size_factor'] = df_sim['size_factors']

        if 'cluster' in set(df_sim.columns):
            adata.obs['cluster'] = df_sim['cluster']        

        # Run test on all of the genes
        print("Running joint test...")
        joint_p_val, t_obs, t_nulls, keep_inds, spotwise_t_obs, spotwise_t_nulls, spot_p_vals = spatialcorr.run_test(
            adata,
            gene_names,
            sigma,
            cond_key=condition_col,
            contrib_thresh=contrib_thresh,
            row_key='row',
            col_key='col',
            verbose=1,
            n_procs=n_procs,
            test_between_conds=between_conds
        ) 
        print(f"joint p-value: {joint_p_val}")
        
        # Run pairwise test
        pair_to_pval = {}
        for g1_i in range(n_genes):
            for g2_i in range(n_genes):
                if g1_i >= g2_i:
                    continue
                print(f"Running test on pair ({g1_i},{g2_i})...")
                pair_p_val, t_obs, t_nulls, keep_inds, spotwise_t_obs, spotwise_t_nulls, spot_p_vals = spatialcorr.run_test(
                    adata,
                    [f'gene_{g1_i}', f'gene_{g2_i}'],
                    sigma,
                    cond_key=condition_col,
                    contrib_thresh=contrib_thresh,
                    row_key='row',
                    col_key='col',
                    verbose=1,
                    n_procs=n_procs,
                    test_between_conds=between_conds
                ) 
                print(f"p-value on pair ({g1_i}, {g2_i}): {pair_p_val}")
                pair_to_pval[f'{g1_i},{g2_i}'] = pair_p_val

        # Store outputs
        dataset_to_results[elem_name]['joint_p_value'] = joint_p_val
        dataset_to_results[elem_name]['pairwise_p_values'] = pair_to_pval

    # Write output
    with open(out_f, 'w') as f:
        json.dump(dataset_to_results, f, indent=True)

if __name__ == '__main__':
    main()
