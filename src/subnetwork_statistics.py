import pandas as pd
import numpy as np
import math
import common
import r_wrappers

params = {
    # "sam_gs_s0": 0.154134905336625  # 3.3  # SAM_GS tuning parameter TODO: chosen to minimize the coefficient of variation
}


def cidrgn1_statistic(
    gene_expr_A,
    gene_expr_B,
    adj_matrix_A,
    adj_matrix_B,
    subnetwork,
    adjusted_regulatory_dissimilarity_norm_factor=1,
):
    """
    Calculate the CIDRGN1 statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)    - adjusted_regulatory_dissimilarity_norm_factor: normalization factor for the adjusted regulatory dissimilarity (default 1)
      (should be computed based on permutation analysis)

    Returns:
    - CIDR-GN1 statistic
    """

    adjusted_regulatory_dissimilarity = adjusted_regulatory_dissimilarity(
        gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B, subnetwork
    )

    sam_gs_stat = sam_gs_statistic(gene_expr_A, gene_expr_B, subnetwork)

    return (
        adjusted_regulatory_dissimilarity
        * adjusted_regulatory_dissimilarity_norm_factor
        + sam_gs_stat
    )


def cidrgn2_statistic(
    gene_expr_A,
    gene_expr_B,
    adj_matrix_A,
    adj_matrix_B,
    subnetwork,
    regulatory_dissimilarity_norm_factor=1,
    topological_dissimilarity_norm_factor=1,
):
    """
    Calculate the CIDRGN2 statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)    - regulatory_dissimilarity_norm_factor: normalization factor for the adjusted regulatory dissimilarity
        (default 1) (should be computed based on permutation analysis)
    - topological_dissimilarity_norm_factor: normalization factor for the topological dissimilarity
        (default 1) (should be computed based on permutation analysis)

    Returns:
    - CIDR-GN1 statistic
    """

    regulatory_dissimilarity = regulatory_dissimilarity(
        gene_expr_A, gene_expr_B, subnetwork
    )

    topological_dissimilarity = topological_dissimilarity(
        adj_matrix_A,
        adj_matrix_B,
        subnetwork,
        common.TopologicalSimilarityMeasure.BASIC,
    )

    sam_gs_stat = sam_gs_statistic(gene_expr_A, gene_expr_B, subnetwork)

    return (
        regulatory_dissimilarity * regulatory_dissimilarity_norm_factor
        + topological_dissimilarity * topological_dissimilarity_norm_factor
        + sam_gs_stat
    )


def regulatory_dissimilarity(
    gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B, subnetwork
):
    """
    Calculate the regulatory dissimilarity statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)
    - subnetwork: list of gene_id for genes in the subnetwork

    Returns:
    - regulatory dissimilarity
    """

    # TODO: should this be symmetric?
    assert (
        adj_matrix_A.shape[0]
        == adj_matrix_A.shape[1]
        == adj_matrix_B.shape[0]
        == adj_matrix_B.shape[1]
    )

    # Compute regulatory dissimilarity per gene
    regulatory_dissimilarity = compute_regulatory_dissimilarities_per_gene(
        gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B
    )

    # Regulatory dissimilarity is averaged across genes in the subnetwork
    total_regulatory_dissimilarity = regulatory_dissimilarity.loc[subnetwork].sum()
    assert total_regulatory_dissimilarity.shape[0] == 1
    return total_regulatory_dissimilarity.iloc[0] / len(subnetwork)


def compute_weighted_regulatory_effect_matrix(adj_matrix, gene_expr):
    """
    Computes the weightedregulatory effect matrix for a given adjacency matrix and gene expression data.

    Parameters:
    - gene_expr: df containing gene expression data (rows are genes, columns are samples)
    - adj_matrix: df containing regression coefficient for the effect of the col gene on row gene for
        (index and column names are gene_id)

    Returns:
    - df containing the weighted regulatory effect matrix
        (the value at (i, j) is the weighted effect of gene j on gene i)
    """

    regulatory_effect_matrix = pd.DataFrame(
        float(0), index=adj_matrix.index, columns=adj_matrix.columns
    )
    for j_gene in regulatory_effect_matrix.columns:
        mean_expr_j_gene = gene_expr.loc[j_gene].mean()
        for i_gene in regulatory_effect_matrix.index:
            regulatory_effect_matrix.loc[i_gene, j_gene] = (
                adj_matrix.loc[i_gene, j_gene] * mean_expr_j_gene
            )
    return regulatory_effect_matrix


def compute_regulatory_dissimilarities_per_gene(
    gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B
):
    """
    Calculate the regulatory dissimilarity of each gene, given the weighted regulatory effect matrices

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)

    Returns:
    - df containing regulatory dissimilarity of each gene (rows are genes, column is "regulatory_dissimilarity")
    """

    # Compute weighted regulatory effect matrices
    regulatory_effect_A = compute_weighted_regulatory_effect_matrix(
        adj_matrix_A, gene_expr_A
    )
    regulatory_effect_B = compute_weighted_regulatory_effect_matrix(
        adj_matrix_B, gene_expr_B
    )

    regulatory_dissimilarity = pd.DataFrame(
        float(0), index=regulatory_effect_A.index, columns=["regulatory_dissimilarity"]
    )
    # Compute regulatory dissimilarity per gene in network
    for gene_id in regulatory_effect_A.index:
        regulatory_effect_vector_A = list(regulatory_effect_A.loc[:, gene_id]) + list(
            regulatory_effect_A.loc[gene_id, :]
        )
        regulatory_effect_vector_B = list(regulatory_effect_B.loc[:, gene_id]) + list(
            regulatory_effect_B.loc[gene_id, :]
        )

        # TODO: This will need some adjustment when we have regulators and targets separated (l2-norm and num_genes should be adjusted)
        num_genes = regulatory_effect_A.shape[0]
        assert (
            len(regulatory_effect_vector_A)
            == len(regulatory_effect_vector_B)
            == 2 * num_genes
        )
        difference_vector = np.array(
            [
                a - b
                for a, b in zip(regulatory_effect_vector_A, regulatory_effect_vector_B)
            ]
        )
        gene_regulatory_dissimilarity = np.sum(difference_vector**2) / len(
            regulatory_effect_vector_A
        )

        regulatory_dissimilarity.loc[gene_id] = gene_regulatory_dissimilarity

    return regulatory_dissimilarity


def topological_dissimilarity(
    adj_matrix_A, adj_matrix_B, subnetwork, topological_dissimilarity_measure
):
    """
    Calculate the topological dissimilarity statistic for a given subnetwork.
    Mathematical definition in README.

    Parameters:
    - adj_matrix_A: df of undirected regulatory network for phenotype A (index and column names are gene_id)
    - adj_matrix_B: df of undirected regtulatory network phenotype B (index and column names are gene_id)
    - subnetwork: list of gene_id for genes in the subnetwork
    - topological_dissimilarity_measure: common.TopologicalDissimilarityMeasure to use

    Returns:
    - topological dissimilarity for the given subnetwork
    """

    # Adjacency matrices should be symmetric
    assert (
        adj_matrix_A.shape[0]
        == adj_matrix_A.shape[1]
        == adj_matrix_B.shape[0]
        == adj_matrix_B.shape[1]
    )
    assert np.allclose(adj_matrix_A, adj_matrix_A.T)
    assert np.allclose(adj_matrix_B, adj_matrix_B.T)

    # Subnetwork genes in adjacency matrices
    for gene in subnetwork:
        assert gene in adj_matrix_A.index, f"Gene {gene} not in adj_matrix_A"
        assert gene in adj_matrix_B.index, f"Gene {gene} not in adj_matrix_B"

    # Compute topolgoical dissimilarity per gene
    if (
        topological_dissimilarity_measure
        == common.TopologicalDissimilarityMeasure.BASIC
    ):
        df_topological = compute_basic_topological_dissimilarity_per_gene(
            adj_matrix_A, adj_matrix_B
        )
    else:
        raise NotImplementedError(
            f"Topological similarity measure {topological_dissimilarity_measure} not implemented"
        )

    # Topolgoical dissimilarity is averaged across genes in the subnetwork
    total_topological_dissimilarity = df_topological.loc[subnetwork].sum()
    assert total_topological_dissimilarity.shape[0] == 1
    return total_topological_dissimilarity.iloc[0] / len(subnetwork)


def compute_basic_topological_dissimilarity_per_gene(adj_matrix_A, adj_matrix_B):
    """
    Calculate the topological dissimilarity statistic for each gene. Based on the basic
    definition in Cidrgn1. Mathematical definition in README.

    Parameters:
    - adj_matrix_A: df of undirected regulatory network for phenotype A (index and column names are gene_id)
    - adj_matrix_B: df of undirected regtulatory network phenotype B (index and column names are gene_id)

    Returns:
    - df containing topological dissimilarity per gene (rows are genes, column is "topological_dissimilarity")
    """

    assert all(adj_matrix_A.index == adj_matrix_B.index)
    assert all(adj_matrix_A.columns == adj_matrix_B.columns)

    topological_dissimilarity = pd.DataFrame(
        float(0), index=adj_matrix_A.index, columns=["topological_dissimilarity"]
    )

    for gene_j in adj_matrix_A.index:
        neighborhood_A = set(
            adj_matrix_A.loc[gene_j][adj_matrix_A.loc[gene_j] == 1].index
        )
        neighborhood_B = set(
            adj_matrix_B.loc[gene_j][adj_matrix_B.loc[gene_j] == 1].index
        )
        neighborhood_similarity = len(
            neighborhood_A.intersection(neighborhood_B)
        ) / len(neighborhood_A.union(neighborhood_B))
        jacard_distance = 1 - neighborhood_similarity
        topological_dissimilarity.loc[gene_j] = jacard_distance

    return topological_dissimilarity


def adjusted_regulatory_dissimilarity(
    gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B, subnetwork
):  # TODO
    """
    Calculate the adjusted regulatory dissimilarity statistic for a given subnetwork.
    Mathematical definition in README. Assuming that a non-zero position in the adjacency matrix
    means that there is an edge between the two genes. (TODO: verify this)

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)
    - subnetwork: list of gene_id for genes in the subnetwork

    Returns:
    - regulatory dissimilarity
    """

    assert (
        adj_matrix_A.shape[0]
        == adj_matrix_A.shape[1]
        == adj_matrix_B.shape[0]
        == adj_matrix_B.shape[1]
    )

    binary_adj_matrix_A = adj_matrix_A.applymap(lambda x: 1 if x != 0 else 0)
    binary_adj_matrix_B = adj_matrix_B.applymap(lambda x: 1 if x != 0 else 0)
    assert np.allclose(binary_adj_matrix_A, adj_matrix_A.T)
    assert np.allclose(binary_adj_matrix_B, binary_adj_matrix_B.T)

    df_topological_dissimilarity = compute_basic_topological_dissimilarity_per_gene(
        binary_adj_matrix_A, binary_adj_matrix_B
    )
    df_regulatory_dissimilarity = compute_regulatory_dissimilarities_per_gene(
        gene_expr_A, gene_expr_B, adj_matrix_A, adj_matrix_B
    )

    df_dissimilarities = df_topological_dissimilarity.merge(df_regulatory_dissimilarity)
    df_dissimilarities["adjusted_regulatory_dissimilarity"] = (
        df_dissimilarities["regulatory_dissimilarity"]
        + df_dissimilarities["topological_dissimilarity"]
        * df_dissimilarities["regulatory_dissimilarity"]
    )

    # Adjusted regulatory dissimilarity is averaged across genes in the subnetwork
    total_adjusted_regulatory_dissimilarity = df_dissimilarities.loc[
        subnetwork, "adjusted_regulatory_dissimilarity"
    ].sum()
    assert total_adjusted_regulatory_dissimilarity.shape[0] == 1
    return total_adjusted_regulatory_dissimilarity.iloc[0] / len(subnetwork)


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


def test_topological_dissimilarity():
    """
    Performs basic sanity checks for the topological dissimilarity implementation
    """
    subnetwork = ["gene_1", "gene_2", "gene_3"]

    adj_matrix_A = pd.read_csv(
        "data/test_data/topological_dissimilarity_adj_matrix_A.tsv",
        sep="\t",
        index_col="gene_id",
    )
    adj_matrix_similar = pd.read_csv(
        "data/test_data/topological_dissimilarity_adj_matrix_similar.tsv",
        sep="\t",
        index_col="gene_id",
    )
    adj_matrix_different = pd.read_csv(
        "data/test_data/topological_dissimilarity_adj_matrix_different.tsv",
        sep="\t",
        index_col="gene_id",
    )

    # Topological dissimilarity between A and A should be 0
    topological_dissimilarity_A = topological_dissimilarity(
        adj_matrix_A,
        adj_matrix_A,
        subnetwork,
        common.TopologicalDissimilarityMeasure.BASIC,
    )
    assert np.allclose(topological_dissimilarity_A, 0)

    # Topological dissimilarity between A and similar should be less than between A and different
    topological_dissimilarity_similar = topological_dissimilarity(
        adj_matrix_A,
        adj_matrix_similar,
        subnetwork,
        common.TopologicalDissimilarityMeasure.BASIC,
    )
    topological_dissimilarity_different = topological_dissimilarity(
        adj_matrix_A,
        adj_matrix_different,
        subnetwork,
        common.TopologicalDissimilarityMeasure.BASIC,
    )
    assert topological_dissimilarity_similar < topological_dissimilarity_different
    assert topological_dissimilarity_similar > 0
    assert topological_dissimilarity_different < 1


def test_regulatory_dissimilarity():
    """
    Performs basic sanity checks for the regulatory dissimilarity implementation
    """
    subnetwork = ["gene_1", "gene_2", "gene_3"]

    adj_matrix_A = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_adj_matrix_A.tsv",
        sep="\t",
        index_col="gene_id",
    )
    adj_matrix_different = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_adj_matrix_different_regulation.tsv",
        sep="\t",
        index_col="gene_id",
    )
    adj_matrix_similar = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_adj_matrix_similar_regulation.tsv",
        sep="\t",
        index_col="gene_id",
    )
    expr_A = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_expr_A.tsv",
        sep="\t",
        index_col="gene_id",
    )
    expr_different = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_expr_different_expr.tsv",
        sep="\t",
        index_col="gene_id",
    )
    expr_similar = pd.read_csv(
        "data/test_data/regulatory_dissimilarity_expr_similar_expr.tsv",
        sep="\t",
        index_col="gene_id",
    )

    # Regulatory dissimilarity between A and A should be 0
    regulatory_dissimilarity_A = regulatory_dissimilarity(
        expr_A,
        expr_A,
        adj_matrix_A,
        adj_matrix_A,
        subnetwork,
    )
    assert np.allclose(regulatory_dissimilarity_A, 0)

    # Regulatory dissimilarity between A and similar regulatory network
    # should be less than between A and different regulatory network
    regulatory_dissimilarity_similar_regulation = regulatory_dissimilarity(
        expr_A,
        expr_A,
        adj_matrix_A,
        adj_matrix_similar,
        subnetwork,
    )
    regulatory_dissimilarity_different_regulation = regulatory_dissimilarity(
        expr_A,
        expr_A,
        adj_matrix_A,
        adj_matrix_different,
        subnetwork,
    )
    assert (
        regulatory_dissimilarity_similar_regulation
        < regulatory_dissimilarity_different_regulation
    )
    assert regulatory_dissimilarity_similar_regulation > 0
    assert (
        regulatory_dissimilarity_different_regulation < 1
    )  # TODO: This probably depends on the values we expect in the regulation adjacency matrix

    # Regulatory dissimilarity between A and similar gene expression
    # should be less than between A and different gene expression
    regulatory_dissimilarity_similar_expr = regulatory_dissimilarity(
        expr_A,
        expr_similar,
        adj_matrix_A,
        adj_matrix_A,
        subnetwork,
    )
    regulatory_dissimilarity_different_expr = regulatory_dissimilarity(
        expr_A,
        expr_different,
        adj_matrix_A,
        adj_matrix_A,
        subnetwork,
    )
    assert (
        regulatory_dissimilarity_similar_expr < regulatory_dissimilarity_different_expr
    )
    assert regulatory_dissimilarity_similar_expr > 0
    assert regulatory_dissimilarity_different_expr < 1


if __name__ == "__main__":
    test_GSCA_implementation()
    test_sam_gs_implementation()
    test_topological_dissimilarity()
    test_regulatory_dissimilarity()
