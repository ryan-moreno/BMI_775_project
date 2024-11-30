import pandas as pd
import numpy as np
import math
import common
import r_wrappers

params = {
    # "sam_gs_s0": 0.154134905336625  # 3.3  # SAM_GS tuning parameter TODO: chosen to minimize the coefficient of variation
}


def gsca_statistic(gene_expr_A, gene_expr_B, subnetwork):
    """
    Calculate the GSCA statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - subnetwork: list of gene_id for genes in the subnetwork

    Returns:
    - GSCA statistic
    """

    num_node_pairs = math.comb(len(subnetwork), 2)

    summed_dispersion_correlations = 0

    for gene in subnetwork:
        for gene_2 in subnetwork:
            if gene == gene_2:
                continue

            # Compute Pearson correlations for gene pair
            correlation_A = gene_expr_A.loc[gene].corr(gene_expr_A.loc[gene_2])
            correlation_B = gene_expr_B.loc[gene].corr(gene_expr_B.loc[gene_2])

            summed_dispersion_correlations += (correlation_A - correlation_B) ** 2

    gsca_statistic = math.sqrt(summed_dispersion_correlations / num_node_pairs)
    return gsca_statistic


def sam_gs_statistic(gene_expr_A, gene_expr_B, subnetwork):
    """
    Calculate the SAM-GS statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - subnetwork: list of gene_id for genes in the subnetwork

    Returns:
    - SAM-GS statistic
    """

    return r_wrappers.sam_gs(gene_expr_A, gene_expr_B, subnetwork)

    sam_gs_total = 0

    for gene in subnetwork:
        gene_mean_expr_A = gene_expr_A.loc[gene].mean()
        gene_mean_expr_B = gene_expr_B.loc[gene].mean()

        # Compute gene-specific scatter
        scaling_factor_numerator = 1 / gene_expr_A.shape[1] + 1 / gene_expr_B.shape[1]
        scaling_factor_denominator = gene_expr_A.shape[1] + gene_expr_B.shape[1] - 2
        scaling_factor = scaling_factor_numerator / scaling_factor_denominator
        unscaled_gene_scatter = 0
        for sample in gene_expr_A.columns:
            unscaled_gene_scatter += (
                gene_expr_A.loc[gene, sample] - gene_mean_expr_A
            ) ** 2
        for sample in gene_expr_B.columns:
            unscaled_gene_scatter += (
                gene_expr_B.loc[gene, sample] - gene_mean_expr_B
            ) ** 2
        gene_scatter = math.sqrt(unscaled_gene_scatter * scaling_factor)

        # Compute gene-specific portion of SAM_GS
        t_like_stat = (gene_mean_expr_A - gene_mean_expr_B) / (
            gene_scatter + params["sam_gs_s0"]
        )
        sam_gs_total += t_like_stat**2

    return sam_gs_total


def test_GSCA_implementation():
    """
    Since GSCA R implementation wasn't working, testing this implementation on their data
    """
    example_data = r_wrappers.load_gsca_example_data()

    calculated_stats = []

    for i, gene_set in enumerate(example_data["gene_sets"]):

        stat = gsca_statistic(
            example_data["gene_expr_A"], example_data["gene_expr_B"], gene_set
        )
        calculated_stats.append(stat)

    # TODO: Currently failing this test. The order is correct, but the values are slightly off
    assert np.allclose(
        calculated_stats, example_data["expected_stats"], rtol=1e-3
    ), f"GSCA statistic implementation incorrect. Expected {example_data['expected_stats']} but got {calculated_stats}"


def test_sam_gs_implementation():
    """
    Tests my implementation against their R code.
    Requires gene expression data to be present for "data/test_data"
    (derived from simulated_gene_expression/sim_gt-cluster_subnetworksize-10_numsamples-50)
    """

    # Load gene expression data
    gene_expr_A_path = "data/test_data/gene_expression_A.tsv"
    gene_expr_B_path = "data/test_data/gene_expression_B.tsv"
    gene_expr_A = pd.read_csv(gene_expr_A_path, sep="\t", index_col="gene_id")
    gene_expr_B = pd.read_csv(gene_expr_B_path, sep="\t", index_col="gene_id")

    # Load metadata
    metadata_path = "data/simulated_gene_expression/sim_gt-cluster_subnetworksize-10_numsamples-50/metadata.tsv"
    metadata = pd.read_csv(metadata_path, sep="\t")

    their_stat_list = []
    our_stat_list = []

    # Get stats for common subnetworks

    for i in range(1, 5):
        subnetwork_id = f"common_subnetwork_{i}"
        subnetwork_genes = metadata[metadata["subnetwork_id"] == subnetwork_id][
            "gene_id"
        ].tolist()

        our_stat = sam_gs_statistic(gene_expr_A, gene_expr_B, subnetwork_genes)
        their_stat = r_wrappers.sam_gs(gene_expr_A, gene_expr_B, subnetwork_genes)
        their_stat_list.append(their_stat)
        our_stat_list.append(our_stat)

    # Get stats for differential subnetworks
    for i in range(1, 7):
        subnetwork_id = f"differential_subnetwork_{i}"
        subnetwork_genes = metadata[metadata["subnetwork_id"] == subnetwork_id][
            "gene_id"
        ].tolist()

        our_stat = sam_gs_statistic(gene_expr_A, gene_expr_B, subnetwork_genes)
        their_stat = r_wrappers.sam_gs(gene_expr_A, gene_expr_B, subnetwork_genes)
        their_stat_list.append(their_stat)
        our_stat_list.append(our_stat)

    assert np.allclose(
        our_stat_list, their_stat_list, rtol=1e-3
    ), f"SAM-GS statistic implementation incorrect. Expected {their_stat_list} but got {our_stat_list}"


if __name__ == "__main__":
    test_GSCA_implementation()
    test_sam_gs_implementation()

    # TODO: Need to test regulatory_dissimilarity
