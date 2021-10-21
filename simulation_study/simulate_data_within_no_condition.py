import scanpy as sc
import ssl
import sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from os.path import join
import spatialcorr
from simulate_spatial import simulate_data
from simulate_spatial import poisson_lognormal
from simulate_spatial import normalize_deviance_residuals as ndr

def main():
    parser = OptionParser()
    parser.add_option("-b", "--bandwidth", help="Kernel bandwidth")
    parser.add_option("-c", "--cov_strength", help="Kernel bandwidth")
    parser.add_option("-o", "--out_dir", help="Output file")
    parser.add_option("-n", "--n_datasets", help="Number of datasets")
    parser.add_option("-s", "--size_factors_type", help="How to compute size factors: 'static_mean', 'permute', 'spotwise'")
    parser.add_option("-g", "--genes", help="Comma-separated list of genes")
    (options, args) = parser.parse_args()

    ssl._create_default_https_context = ssl._create_unverified_context
    
    # Get parameters
    sigma = float(options.bandwidth)
    cov_strength = float(options.cov_strength)
    n_datasets = int(options.n_datasets)
    gene_1, gene_2 = options.genes.split(",")
    size_factors_type = options.size_factors_type
    out_dir = options.out_dir

    # Load data
    print("Loading dataset...")
    adata_raw = sc.datasets.visium_sge(sample_id="V1_Breast_Cancer_Block_A_Section_1")
    sc.pp.filter_cells(adata_raw, min_counts=5000)
    sc.pp.filter_cells(adata_raw, max_counts=35000)
    adata_raw.X = np.array(adata_raw.X.todense())

    # Compute kernel matrix
    print("Computing kernel matrix...")
    kernel_matrix = spatialcorr.statistical_test._compute_kernel_matrix(
        adata_raw.obs, 
        sigma=sigma, 
        y_col='array_row', 
        x_col='array_col'
    )

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

    # Infer latent parameters from real data
    print(f"Inferring simulation parameters from real data for gene {gene_1}...")
    means_1, vars_1 = poisson_lognormal.fit(
        adata_raw.obs_vector(gene_1),
        size_factors.T[0]
    )
    mean_g1 = np.mean(means_1.squeeze())
    var_g1 = np.mean(vars_1.squeeze())**2
    print(f"Inferring simulation parameters from real data for gene {gene_2}...")
    means_2, vars_2 = poisson_lognormal.fit(
        adata_raw.obs_vector(gene_2),
        size_factors.T[1]
    )
    mean_g2 = np.mean(means_2.squeeze())
    var_g2 = np.mean(vars_2.squeeze())**2

    for dataset_i in range(n_datasets):
        print(f"Generating dataset {dataset_i+1}...")
        # Simulate dataset
        sim_corrs, sim_covs, adata_sim = simulate_data.simulate_pairwise_from_dataset(
            adata_raw,
            corr_mean, 
            gene_1,
            gene_2,
            row_key='array_row',
            col_key='array_col',
            sigma=sigma,
            cov_strength=cov_strength,
            poisson=True,
            size_factors=size_factors,
            mean_g1=mean_g1,
            mean_g2=mean_g2,
            var_g1=var_g1,
            var_g2=var_g2
        )

        # Normalize simulated samples
        sim_X = ndr.normalize(
            np.array(
                adata_sim.X, 
                dtype=np.float64
            ), 
            umi_counts=np.array(
                size_factors.T[0], 
                dtype=np.float64
            )
        )

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
                "size_factors": size_factors.T[0]
            },
            index=adata_raw.obs.index
        )
        df_sim.to_csv(join(out_dir, f"{dataset_i}.tsv"), sep='\t')


if __name__ == '__main__':
    main()
