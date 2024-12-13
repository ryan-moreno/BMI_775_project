import pandas as pd
import numpy as np
import os
import common
import subnetwork_statistics
from sklearn.preprocessing import StandardScaler
import ast

params = {
    "n_permutations": 1000,
}

standardScaler = StandardScaler()


def try_node2vec_parameters(
    adj_matrix_A,
    gene_expression_A,
    adj_matrix_B,
    gene_expression_B,
    subnetworks_of_interest,
    num_permutations,
    existing_perm_test_stats_path=None,
):
    """
    Performs a permutation test to benchmark the performance of the node2vec based statistics.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)
    - subnetworks_of_interest: dictionary with key subnetwork_id and value list of gene_id for genes in the subnetwork
    - num_permutations: number of permutations to perform
    - existing_perm_test_stats_path: path to existing permutation test stats to use as a starting point (and to use those permuted subnetworks)

    Returns:
    - A tuple containing:
        - permutation_test_stats: df containing the test statistics for each subnetwork (indices are subnetwork_id and columns are statistics)
        - p_values: df containing the p-values for each subnetwork (indices are subnetwork_id for subnetworks of interest and
            columns are p-values associated with statistics)
    """

    p_q_node2vec_params = [
        (2, 2, 10)
        # (1, 1, 5),
        # (1, 1, 10),
        # (1, 1, 15),
        # (1, 2, 10),
        # (1, 0.5, 10),
        # (2, 1, 10),
        # (2, 2, 10),
        # (2, 0.5, 10),
        # (0.5, 1, 10),
        # (0.5, 2, 10),
        # (0.5, 0.5, 10),
    ]

    assert (
        adj_matrix_A.index.equals(adj_matrix_A.columns)
        and adj_matrix_A.index.equals(adj_matrix_B.index)
        and adj_matrix_B.index.equals(adj_matrix_B.columns)
        and adj_matrix_A.index.equals(gene_expression_A.index)
        and gene_expression_A.index.equals(gene_expression_B.index)
    ), "Index and column dimensions do not align correctly."

    full_gene_list = list(adj_matrix_A.index)

    # All subnetworks should have the same length
    num_genes = len(subnetworks_of_interest[list(subnetworks_of_interest.keys())[0]])
    assert all(
        len(subnetworks_of_interest[list(subnetworks_of_interest.keys())[i]])
        == num_genes
        for i in range(1, len(subnetworks_of_interest))
    ), "Subnetworks do not have the same size."

    if existing_perm_test_stats_path is None:

        # Generate permuted subnetworks
        permuted_subnetworks = []
        for i in range(num_permutations):
            random_subnetwork = np.random.choice(
                full_gene_list, num_genes, replace=False
            )
            permuted_subnetworks.append(list(random_subnetwork))

        # Initialize dataframe to store test statistics
        # Put subnetworks of interest at the top, followed by permutated subnetworks
        subnetwork_ids = list(subnetworks_of_interest.keys()) + [
            f"permuted_{i}" for i in range(num_permutations)
        ]
        permutation_test_stats = pd.DataFrame(
            index=subnetwork_ids,
            columns=[
                "subnetwork_of_interest",
                "subnetwork",
                "regulatory_dissimilarity",
                "sam_gs",
            ]
            + [
                f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
                for (p, q, n_neighbors) in p_q_node2vec_params
            ]
            + [
                f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
                for (p, q, n_neighbors) in p_q_node2vec_params
            ],
        )

        permutation_test_stats["subnetwork_of_interest"] = [True] * len(
            subnetworks_of_interest.keys()
        ) + [False] * len(permuted_subnetworks)
        permutation_test_stats["subnetwork"] = (
            list(subnetworks_of_interest.values()) + permuted_subnetworks
        )
    else:
        permutation_test_stats = pd.read_csv(
            existing_perm_test_stats_path,
            sep="\t",
            index_col=0,
            converters={"subnetwork": ast.literal_eval},
        )
        for p, q, n_neighbors in p_q_node2vec_params:
            permutation_test_stats[
                f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            ] = np.nan
            permutation_test_stats[
                f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            ] = np.nan

    # Compute node2vec embeddings for both graphs
    node2vec_embeddings_A = []
    node2vec_embeddings_B = []
    for p, q, n_neighbors in p_q_node2vec_params:
        embeddings_A, embeddings_B = subnetwork_statistics.compute_node2vec_embeddings(
            adj_matrix_A, adj_matrix_B, p, q
        )
        node2vec_embeddings_A.append(embeddings_A)
        node2vec_embeddings_B.append(embeddings_B)

    # Calculate test statistics for each subnetwork (subnetwork of interest & permuted ones)
    for subnetwork_id in permutation_test_stats.index:
        print(f"Calculating test statistics for subnetwork {subnetwork_id}...")

        subnetwork = permutation_test_stats.loc[subnetwork_id]["subnetwork"]

        # When reading in from existing permutation test stats, some of these were formatted incorrectly
        if len(subnetwork) == 1 and not "," in subnetwork[0]:
            subnetwork = subnetwork[0].split("gene_")[1:]
            subnetwork = [f"gene_{gene}" for gene in subnetwork]

        assert len(subnetwork) == num_genes, "Subnetworks do not have the same size."

        is_subnetwork_of_interest = permutation_test_stats.loc[subnetwork_id][
            "subnetwork_of_interest"
        ]
        regulatory_dissimilarity = subnetwork_statistics.regulatory_dissimilarity(
            gene_expression_A, gene_expression_B, adj_matrix_A, adj_matrix_B, subnetwork
        )

        topological_dissimilarities_node2vec = []
        adjusted_regulatory_dissimilarities_node2vec = []
        for i, (p, q, n_neighbors) in enumerate(p_q_node2vec_params):
            topological_dissimilarity_node2vec = (
                subnetwork_statistics.topological_dissimilarity(
                    adj_matrix_A,
                    adj_matrix_B,
                    subnetwork,
                    common.TopologicalDissimilarityMeasure.NODE2VEC,
                    node2vec_embeddings_A[i],
                    node2vec_embeddings_B[i],
                    n_neighbors,
                )
            )
            topological_dissimilarities_node2vec.append(
                topological_dissimilarity_node2vec
            )

            adjusted_regulatory_dissimilarity_node2vec = (
                subnetwork_statistics.adjusted_regulatory_dissimilarity(
                    gene_expression_A,
                    gene_expression_B,
                    adj_matrix_A,
                    adj_matrix_B,
                    common.TopologicalDissimilarityMeasure.NODE2VEC,
                    subnetwork,
                    node2vec_embeddings_A[i],
                    node2vec_embeddings_B[i],
                    n_neighbors,
                )
            )
            adjusted_regulatory_dissimilarities_node2vec.append(
                adjusted_regulatory_dissimilarity_node2vec
            )
        sam_gs = subnetwork_statistics.sam_gs_statistic(
            gene_expression_A, gene_expression_B, subnetwork
        )

        subnetwork_stat_dictionary = {
            "subnetwork_of_interest": is_subnetwork_of_interest,
            "subnetwork": subnetwork,
            "regulatory_dissimilarity": regulatory_dissimilarity,
            "sam_gs": sam_gs,
        }

        for i, (p, q, n_neighbors) in enumerate(p_q_node2vec_params):
            subnetwork_stat_dictionary[
                f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            ] = topological_dissimilarities_node2vec[i]
            subnetwork_stat_dictionary[
                f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            ] = adjusted_regulatory_dissimilarities_node2vec[i]

        permutation_test_stats.loc[subnetwork_id] = subnetwork_stat_dictionary

    # Normalize cidrgn related statistics
    permutation_test_stats["regulatory_dissimilarity_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["regulatory_dissimilarity"]]
        )
    )
    for i, (p, q, n_neighbors) in enumerate(p_q_node2vec_params):
        permutation_test_stats[
            f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}_norm"
        ] = standardScaler.fit_transform(
            permutation_test_stats[
                [f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"]
            ]
        )
        permutation_test_stats[
            f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}_norm"
        ] = standardScaler.fit_transform(
            permutation_test_stats[
                [
                    f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
                ]
            ]
        )

    # Calculate cidrgn
    for i, (p, q, n_neighbors) in enumerate(p_q_node2vec_params):
        permutation_test_stats[f"cidrgn1_node2vec_p_{p}_q_{q}_m_{n_neighbors}"] = (
            permutation_test_stats[
                f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}_norm"
            ]
            + permutation_test_stats["sam_gs"]
        )
        permutation_test_stats[f"cidrgn2_node2vec_p_{p}_q_{q}_m_{n_neighbors}"] = (
            permutation_test_stats["regulatory_dissimilarity_norm"]
            + permutation_test_stats[
                f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}_norm"
            ]
            + permutation_test_stats["sam_gs"]
        )

    # Compute p-values
    p_values = pd.DataFrame(
        index=list(subnetworks_of_interest.keys()),
        columns=[
            "sam_gs",
            "regulatory_dissimilarity",
        ]
        + [
            f"cidrgn1_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            for (p, q, n_neighbors) in p_q_node2vec_params
        ]
        + [
            f"cidrgn2_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            for (p, q, n_neighbors) in p_q_node2vec_params
        ]
        + [
            f"topological_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            for (p, q, n_neighbors) in p_q_node2vec_params
        ]
        + [
            f"adjusted_regulatory_dissimilarity_node2vec_p_{p}_q_{q}_m_{n_neighbors}"
            for (p, q, n_neighbors) in p_q_node2vec_params
        ],
    )

    for subnetwork_id in p_values.index:
        for statistic_type in p_values.columns:
            val_of_interest = permutation_test_stats.loc[subnetwork_id][statistic_type]
            num_random_trials_more_extreme = sum(
                val_of_interest
                <= permutation_test_stats[
                    permutation_test_stats["subnetwork_of_interest"] == False
                ][statistic_type]
            )

            p_values.loc[subnetwork_id][statistic_type] = (
                num_random_trials_more_extreme / num_permutations
            )

    return permutation_test_stats, p_values


def benchmark_permutation_test(
    adj_matrix_A,
    gene_expression_A,
    adj_matrix_B,
    gene_expression_B,
    subnetworks_of_interest,
    num_permutations,
):
    """
    Performs a permutation test to benchmark the performance of a given statistic method.

    Parameters:
    - gene_expr_A: df containing gene expression data for phenotype A (rows are genes, columns are samples)
    - gene_expr_B: df containing gene expression data for phenotype B (rows are genes, columns are samples)
    - adj_matrix_A: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype A (index and column names are gene_id)
    - adj_matrix_B: df containing regression coefficient for the effect of the col gene on row gene for
        phenotype B (index and column names are gene_id)
    - subnetworks_of_interest: dictionary with key subnetwork_id and value list of gene_id for genes in the subnetwork
    - num_permutations: number of permutations to perform

    Returns:
    - A tuple containing:
        - permutation_test_stats: df containing the test statistics for each subnetwork (indices are subnetwork_id and columns are statistics)
        - p_values: df containing the p-values for each subnetwork (indices are subnetwork_id for subnetworks of interest and
            columns are p-values associated with statistics)
    """

    assert (
        adj_matrix_A.index.equals(adj_matrix_A.columns)
        and adj_matrix_A.index.equals(adj_matrix_B.index)
        and adj_matrix_B.index.equals(adj_matrix_B.columns)
        and adj_matrix_A.index.equals(gene_expression_A.index)
        and gene_expression_A.index.equals(gene_expression_B.index)
    ), "Index and column dimensions do not align correctly."

    full_gene_list = list(adj_matrix_A.index)

    # All subnetworks should have the same length
    num_genes = len(subnetworks_of_interest[list(subnetworks_of_interest.keys())[0]])
    assert all(
        len(subnetworks_of_interest[list(subnetworks_of_interest.keys())[i]])
        == num_genes
        for i in range(1, len(subnetworks_of_interest))
    ), "Subnetworks do not have the same size."

    # Generate permuted subnetworks
    permuted_subnetworks = []
    for i in range(num_permutations):
        random_subnetwork = np.random.choice(full_gene_list, num_genes, replace=False)
        permuted_subnetworks.append(random_subnetwork)

    # Initialize dataframe to store test statistics
    # Put subnetworks of interest at the top, followed by permutated subnetworks
    subnetwork_ids = list(subnetworks_of_interest.keys()) + [
        f"permuted_{i}" for i in range(num_permutations)
    ]
    permutation_test_stats = pd.DataFrame(
        index=subnetwork_ids,
        columns=[
            "subnetwork_of_interest",
            "subnetwork",
            "regulatory_dissimilarity",
            "topological_dissimilarity_basic",
            "topological_dissimilarity_node2vec",
            "adjusted_regulatory_dissimilarity_basic",
            "adjusted_regulatory_dissimilarity_node2vec",
            "sam_gs",
            "gsca",
        ],
    )

    permutation_test_stats["subnetwork_of_interest"] = [True] * len(
        subnetworks_of_interest.keys()
    ) + [False] * len(permuted_subnetworks)
    permutation_test_stats["subnetwork"] = (
        list(subnetworks_of_interest.values()) + permuted_subnetworks
    )

    # Compute node2vec embeddings for both graphs
    embeddings_A, embeddings_B = subnetwork_statistics.compute_node2vec_embeddings(
        adj_matrix_A, adj_matrix_B
    )

    # Calculate test statistics for each subnetwork (subnetwork of interest & permuted ones)
    for subnetwork_id in permutation_test_stats.index:
        print(f"Calculating test statistics for subnetwork {subnetwork_id}...")

        subnetwork = permutation_test_stats.loc[subnetwork_id]["subnetwork"]
        assert len(subnetwork) == num_genes, "Subnetworks do not have the same size."

        is_subnetwork_of_interest = permutation_test_stats.loc[subnetwork_id][
            "subnetwork_of_interest"
        ]
        regulatory_dissimilarity = subnetwork_statistics.regulatory_dissimilarity(
            gene_expression_A, gene_expression_B, adj_matrix_A, adj_matrix_B, subnetwork
        )
        topological_dissimilarity_basic = (
            subnetwork_statistics.topological_dissimilarity(
                adj_matrix_A,
                adj_matrix_B,
                subnetwork,
                common.TopologicalDissimilarityMeasure.BASIC,
            )
        )
        topological_dissimilarity_node2vec = (
            subnetwork_statistics.topological_dissimilarity(
                adj_matrix_A,
                adj_matrix_B,
                subnetwork,
                common.TopologicalDissimilarityMeasure.NODE2VEC,
                embeddings_A,
                embeddings_B,
            )
        )
        adjusted_regulatory_dissimilarity_basic = (
            subnetwork_statistics.adjusted_regulatory_dissimilarity(
                gene_expression_A,
                gene_expression_B,
                adj_matrix_A,
                adj_matrix_B,
                common.TopologicalDissimilarityMeasure.BASIC,
                subnetwork,
            )
        )
        adjusted_regulatory_dissimilarity_node2vec = (
            subnetwork_statistics.adjusted_regulatory_dissimilarity(
                gene_expression_A,
                gene_expression_B,
                adj_matrix_A,
                adj_matrix_B,
                common.TopologicalDissimilarityMeasure.NODE2VEC,
                subnetwork,
                embeddings_A,
                embeddings_B,
            )
        )
        sam_gs = subnetwork_statistics.sam_gs_statistic(
            gene_expression_A, gene_expression_B, subnetwork
        )
        gsca = subnetwork_statistics.gsca_statistic(
            gene_expression_A, gene_expression_B, subnetwork
        )

        permutation_test_stats.loc[subnetwork_id] = {
            "subnetwork_of_interest": is_subnetwork_of_interest,
            "subnetwork": subnetwork,
            "regulatory_dissimilarity": regulatory_dissimilarity,
            "topological_dissimilarity_basic": topological_dissimilarity_basic,
            "topological_dissimilarity_node2vec": topological_dissimilarity_node2vec,
            "adjusted_regulatory_dissimilarity_basic": adjusted_regulatory_dissimilarity_basic,
            "adjusted_regulatory_dissimilarity_node2vec": adjusted_regulatory_dissimilarity_node2vec,
            "sam_gs": sam_gs,
            "gsca": gsca,
        }

    # Normalize cidrgn related statistics
    # TODO: Should the normalization be done with all of the subnetworks of interest at once?
    permutation_test_stats["regulatory_dissimilarity_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["regulatory_dissimilarity"]]
        )
    )
    permutation_test_stats["topological_dissimilarity_basic_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["topological_dissimilarity_basic"]]
        )
    )
    permutation_test_stats["topological_dissimilarity_node2vec_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["topological_dissimilarity_node2vec"]]
        )
    )
    permutation_test_stats["adjusted_regulatory_dissimilarity_basic_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["adjusted_regulatory_dissimilarity_basic"]]
        )
    )
    permutation_test_stats["adjusted_regulatory_dissimilarity_node2vec_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["adjusted_regulatory_dissimilarity_node2vec"]]
        )
    )

    # Calculate cidrgn
    permutation_test_stats["cidrgn1"] = (
        permutation_test_stats["adjusted_regulatory_dissimilarity_basic_norm"]
        + permutation_test_stats["sam_gs"]
    )
    permutation_test_stats["cidrgn2"] = (
        permutation_test_stats["regulatory_dissimilarity_norm"]
        + permutation_test_stats["topological_dissimilarity_basic_norm"]
        + permutation_test_stats["sam_gs"]
    )
    permutation_test_stats["cidrgn1_node2vec"] = (
        permutation_test_stats["adjusted_regulatory_dissimilarity_node2vec_norm"]
        + permutation_test_stats["sam_gs"]
    )
    permutation_test_stats["cidrgn2_node2vec"] = (
        permutation_test_stats["regulatory_dissimilarity_norm"]
        + permutation_test_stats["topological_dissimilarity_node2vec_norm"]
        + permutation_test_stats["sam_gs"]
    )

    # Compute p-values
    p_values = pd.DataFrame(
        index=list(subnetworks_of_interest.keys()),
        columns=[
            "sam_gs",
            "gsca",
            "cidrgn1",
            "cidrgn2",
            "cidrgn1_node2vec",
            "cidrgn2_node2vec",
            "regulatory_dissimilarity",
            "topological_dissimilarity_basic",
            "topological_dissimilarity_node2vec",
            "adjusted_regulatory_dissimilarity_basic",
            "adjusted_regulatory_dissimilarity_node2vec",
        ],
    )

    for subnetwork_id in p_values.index:
        for statistic_type in p_values.columns:
            val_of_interest = permutation_test_stats.loc[subnetwork_id][statistic_type]
            num_random_trials_more_extreme = sum(
                val_of_interest
                <= permutation_test_stats[
                    permutation_test_stats["subnetwork_of_interest"] == False
                ][statistic_type]
            )

            p_values.loc[subnetwork_id][statistic_type] = (
                num_random_trials_more_extreme / num_permutations
            )

    return permutation_test_stats, p_values


def permutation_tests_simulated_data():
    """
    Runs a permutation test for each of the simulated data sets.
    Requires the simulated gene expression data and the inferred regulatory networks to be present.
    """
    for graph_type in common.SubnetworkType:
        if graph_type == common.SubnetworkType.RANDOM:
            continue  # The way I have this set up right now, RANDOM and CLUSTER are the same
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"
            )

            # Read in gene expression data
            gene_expr_A = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/gene_expression_A.tsv",
                sep="\t",
                index_col="gene_id",
            )
            gene_expr_B = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/gene_expression_B.tsv",
                sep="\t",
                index_col="gene_id",
            )
            metadata = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/metadata.tsv",
                sep="\t",
                index_col="gene_id",
            )

            # Read in regulatory network
            adj_matrix_A = pd.read_csv(
                f"data/simulated_inferred_networks/{simulation_name}/adj_matrix_A.tsv",
                sep="\t",
                index_col="gene_id",
            )
            adj_matrix_B = pd.read_csv(
                f"data/simulated_inferred_networks/{simulation_name}/adj_matrix_B.tsv",
                sep="\t",
                index_col="gene_id",
            )

            # Get the subnetworks of interest
            subnetwork_ids = metadata.loc[
                metadata["subnetwork_id"] != "remaining_genes", "subnetwork_id"
            ].unique()

            subnetworks_of_interest = {}
            for subnetwork_id in subnetwork_ids:
                subnetworks_of_interest[subnetwork_id] = list(
                    metadata[metadata["subnetwork_id"] == subnetwork_id].index
                )

            # Perform permutation test
            permutation_test_stats, p_values = benchmark_permutation_test(
                adj_matrix_A,
                gene_expr_A,
                adj_matrix_B,
                gene_expr_B,
                subnetworks_of_interest,
                params["n_permutations"],
            )

            # Save results
            folder = f"data/simulated_data_results/{simulation_name}"
            if not os.path.exists(folder):
                os.makedirs(folder)

            permutation_test_stats.to_csv(
                f"{folder}/permutation_test_stats.tsv",
                sep="\t",
            )
            p_values.to_csv(
                f"{folder}/p_values.tsv",
                sep="\t",
            )


def permutation_tests_node2vec_try_pq_params(use_rep=None):
    """
    Runs a permutation test for each of the simulated data sets.
    Requires the simulated gene expression data and the inferred regulatory networks to be present.
    Parameters:
    - use_rep: If not None, will use the permutation stats from this rep as the permutations to use
    """
    for graph_type in common.SubnetworkType:
        if graph_type == common.SubnetworkType.RANDOM:
            continue  # The way I have this set up right now, RANDOM and CLUSTER are the same
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"
            )

            # Read in gene expression data
            gene_expr_A = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/gene_expression_A.tsv",
                sep="\t",
                index_col="gene_id",
            )
            gene_expr_B = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/gene_expression_B.tsv",
                sep="\t",
                index_col="gene_id",
            )
            metadata = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/metadata.tsv",
                sep="\t",
                index_col="gene_id",
            )

            # Read in regulatory network
            adj_matrix_A = pd.read_csv(
                f"data/simulated_inferred_networks/{simulation_name}/adj_matrix_A.tsv",
                sep="\t",
                index_col="gene_id",
            )
            adj_matrix_B = pd.read_csv(
                f"data/simulated_inferred_networks/{simulation_name}/adj_matrix_B.tsv",
                sep="\t",
                index_col="gene_id",
            )

            # Get the subnetworks of interest
            subnetwork_ids = metadata.loc[
                metadata["subnetwork_id"] != "remaining_genes", "subnetwork_id"
            ].unique()

            subnetworks_of_interest = {}
            for subnetwork_id in subnetwork_ids:
                subnetworks_of_interest[subnetwork_id] = list(
                    metadata[metadata["subnetwork_id"] == subnetwork_id].index
                )

            if use_rep is not None:
                existing_perm_test_stats_path = f"data/simulated_data_results/{use_rep}/{simulation_name}/permutation_test_stats.tsv"
            else:
                existing_perm_test_stats_path = None

            # Perform permutation test
            permutation_test_stats, p_values = try_node2vec_parameters(
                adj_matrix_A,
                gene_expr_A,
                adj_matrix_B,
                gene_expr_B,
                subnetworks_of_interest,
                params["n_permutations"],
                existing_perm_test_stats_path=existing_perm_test_stats_path,
            )

            # Save results
            folder = f"data/simulated_data_results/{simulation_name}"
            if not os.path.exists(folder):
                os.makedirs(folder)

            permutation_test_stats.to_csv(
                f"{folder}/node2vec_param_trials_permutation_test_stats.tsv",
                sep="\t",
            )
            p_values.to_csv(
                f"{folder}/node2vec_param_trials_p_values.tsv",
                sep="\t",
            )


if __name__ == "__main__":
    # permutation_tests_simulated_data()
    permutation_tests_node2vec_try_pq_params("rep7")
