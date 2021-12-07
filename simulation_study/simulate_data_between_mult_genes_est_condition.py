import scanpy as sc
import ssl
import sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import pickle
from os.path import join
import json
from collections import defaultdict

import spatialcorr
from simulate_spatial import simulate_data_condition as simulate_data
from simulate_spatial import poisson_lognormal
from simulate_spatial import normalize_deviance_residuals as ndr

def main():
    parser = OptionParser()
    parser.add_option("-o", "--out_dir", help="Output file")
    parser.add_option("-n", "--n_datasets", help="Number of datasets")
    parser.add_option("-s", "--size_factors_type", help="How to compute size factors: 'static_mean', 'permute', 'spotwise'")
    parser.add_option("-i", "--config", help="Simulation config file")
    (options, args) = parser.parse_args()

    ssl._create_default_https_context = ssl._create_unverified_context
    
    # Get parameters
    n_datasets = int(options.n_datasets)
    size_factors_type = options.size_factors_type
    config_f = options.config
    out_dir = options.out_dir

    with open(config_f, 'r') as f:
        config = json.load(f)

    # Load data
    print("Loading dataset...")
    adata_raw = sc.datasets.visium_sge(sample_id="V1_Breast_Cancer_Block_A_Section_1")
    sc.pp.filter_cells(adata_raw, min_counts=5000)
    sc.pp.filter_cells(adata_raw, max_counts=35000)
    adata_raw.X = np.array(adata_raw.X.todense())

    # Add cluster info from config
    assert frozenset(adata_raw.obs.index) == frozenset(config['spot_to_cluster'].keys())
    adata_raw.obs['cluster'] = [
        config['spot_to_cluster'][spot]
        for spot in adata_raw.obs.index
    ]

    # Get genes from config
    the_genes = config['genes']

    # Compute the size-factors
    umi_counts = np.sum(adata_raw.X, axis=1)
    #umi_counts = np.sum(np.array(adata_raw.X.todense()), axis=1)

    if size_factors_type == 'spotwise':
        size_factors = np.full(
            (len(the_genes), adata_raw.X.shape[0]),
            umi_counts
        ).T
    elif size_factors_type == 'permute':
        size_factors = np.full(
            (len(the_genes), adata_raw.X.shape[0]),
            umi_counts
        ).T
        size_factors = np.random.permutation(size_factors)
    elif size_factors_type == 'static_mean':
        size_factors = np.full(
            (len(the_genes), adata_raw.X.shape[0]),
            np.mean(umi_counts).astype(int)
        ).T

    # Map each cluster to its indices
    clust_to_inds = defaultdict(lambda: [])
    for ind, ct in enumerate(adata_raw.obs['cluster']):
        clust_to_inds[ct].append(ind)

    # Estimate latent mean and variance within each cluster
    gene_to_clust_to_mean = defaultdict(lambda: {})
    gene_to_clust_to_var = defaultdict(lambda: {})
    for gene in set(the_genes):
        for clust in sorted(set(adata_raw.obs['cluster'])):
            print(f"Inferring latent mean and variance for gene {gene} in cluster {clust}...")
            inds = clust_to_inds[clust]
            means, varss = poisson_lognormal.fit(
                adata_raw.obs_vector(gene)[inds],
                size_factors.T[0][inds]
            )
            mean = np.mean(means.squeeze())
            varr = np.mean(varss.squeeze())**2
            gene_to_clust_to_mean[gene][clust] = mean
            gene_to_clust_to_var[gene][clust] = varr

    # Create spotwise covariance matrices
    spotwise_covs = []
    for s_i, spot in enumerate(adata_raw.obs.index):
        clust = config['spot_to_cluster'][spot]
        spot_corr = config['cluster_to_correlation'][clust]
        cov_mat = np.zeros((len(the_genes), len(the_genes)))
        for g1_i, g1 in enumerate(the_genes):
            for g2_i, g2 in enumerate(the_genes):
                if g1_i == g2_i:
                    varr = gene_to_clust_to_var[g1][clust]
                    cov_mat[g1_i][g2_i] = varr
                else:
                    varr_1 = gene_to_clust_to_var[g1][clust]
                    varr_2 = gene_to_clust_to_var[g2][clust]
                    cov = np.sqrt(varr_1 * varr_2) * spot_corr
                    cov_mat[g1_i][g2_i] = cov
        spotwise_covs.append(cov_mat)

    # Simulate a collection of datasets
    for dataset_i in range(n_datasets):
        print(f"Generating dataset {dataset_i+1}...")
        # Simulate dataset
        sim_corrs, sim_covs, adata_sim = simulate_data.simulate_gene_set_from_dataset_precomputed_cov( 
            adata_raw,
            the_genes,
            spotwise_covs,
            row_key='array_row',
            col_key='array_col',
            clust_key='cluster',
            poisson=True,
            size_factors=size_factors,
            gene_to_clust_to_mean=gene_to_clust_to_mean,
            gene_to_clust_to_var=gene_to_clust_to_var
        )

        print(adata_sim.X)
        
        # Normalize simulated samples
        sim_X = ndr.normalize(
            np.array(adata_sim.X, dtype=np.float64), 
            umi_counts=np.array(size_factors.T[0], dtype=np.float64)
        )

        # Write dataset
        data={
            "row": adata_raw.obs['array_row'],
            "col": adata_raw.obs['array_col'],
            "correlation": [str(x.tolist()) for x in sim_corrs],
            "covariance": [str(x.tolist()) for x in sim_covs],
            "size_factors": size_factors.T[0],
            "cluster": list(adata_raw.obs['cluster'])
        }
        for gene_i in range(adata_sim.X.shape[1]):
            data[f'count_gene_{gene_i}'] = adata_sim.X.T[gene_i]
        for gene_i in range(sim_X.shape[1]):
            data[f'dev_res_gene_{gene_i}'] = sim_X.T[gene_i]
        df_sim = pd.DataFrame(
            data=data,
            index=adata_raw.obs.index
        )
        df_sim.to_csv(join(out_dir, f"{dataset_i}.tsv"), sep='\t')


if __name__ == '__main__':
    main()
