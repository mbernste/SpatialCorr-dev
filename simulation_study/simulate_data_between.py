import scanpy as sc
import ssl
import sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import json
from os.path import join

from simulate_spatial import simulate_data
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

    # Set the mean correlation across the slide to zero
    corr_mean = 0

    # Compute the size-factors
    umi_counts = np.sum(adata_raw.X, axis=1)

    if size_factors_type == 'spotwise':
        size_factors = np.full(
            (2, adata_raw.X.shape[0]),
            umi_counts
        ).T
    elif size_factors_type == 'permute':
        size_factors = np.full(
            (2, adata_raw.X.shape[0]),
            umi_counts
        ).T
        size_factors = np.random.permutation(size_factors)
    elif size_factors_type == 'static_mean':
        size_factors = np.full(
            (2, adata_raw.X.shape[0]),
            np.mean(umi_counts)
        ).T


    all_genes = set(config['cluster_to_gene_stat_1'].values()) | set(config['cluster_to_gene_stat_2'].values())
    gene_to_mean = {}
    gene_to_var = {}
    for gene in all_genes:
        means, varss = poisson_lognormal.fit(
            adata_raw.obs_vector(gene),
            size_factors.T[0]
        )
        mean = np.mean(means.squeeze())
        varr = np.mean(varss.squeeze())**2
        gene_to_mean[gene] = mean
        gene_to_var[gene] = varr

    for dataset_i in range(n_datasets):
        print(f"Generating dataset {dataset_i+1}...")
        
        # Simulate dataset
        sim_corrs, sim_covs, adata_sim =simulate_data.simulate_pairwise_between_groups_from_dataset(
            adata_raw,
            clust_to_corr=config['cluster_to_correlation'],
            clust_to_gene_stat_1=config['cluster_to_gene_stat_1'],
            clust_to_gene_stat_2=config['cluster_to_gene_stat_2'],
            clust_key='cluster',
            row_key='array_row',
            col_key='array_col',
            poisson=True,
            size_factors=size_factors,
            gene_to_mean=gene_to_mean,
            gene_to_var=gene_to_var
        )

        # Normalize simulated samples
        sim_X = ndr.normalize(adata_sim.X, umi_counts=np.array(size_factors.T[0], dtype=np.float64))

        # Write dataset
        df_sim = pd.DataFrame(
            data={
                "row": adata_raw.obs['array_row'],
                "col": adata_raw.obs['array_col'],
                "count_gene_1": adata_sim.X.T[0],
                "count_gene_2": adata_sim.X.T[1],
                "dev_res_gene_1": sim_X.T[0],
                "dev_res_gene_2": sim_X.T[1],
                "correlation": sim_corrs,
                "covariance": sim_covs,
                "size_factors": size_factors.T[0],
                "cluster": list(adata_raw.obs['cluster'])
            },
            index=adata_raw.obs.index
        )
        df_sim.to_csv(join(out_dir, f"{dataset_i}.tsv"), sep='\t')

if __name__ == '__main__':
    main()
