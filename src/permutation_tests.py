import pandas as pd
import numpy as np
import os
import common
import subnetwork_statistics
from sklearn.preprocessing import StandardScaler

params = {
    "n_permutations": 1000,
}

standardScaler = StandardScaler()


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
    - A dictionary containing the following:
        - p_value: The p-value of the permutation test
        - test_statistic: The test statistic of the permutation test
        - null_distribution: The null distribution of the permutation test
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
            "topological_dissimilarity",
            "adjusted_regulatory_dissimilarity",
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
        topological_dissimilarity = subnetwork_statistics.topological_dissimilarity(
            adj_matrix_A,
            adj_matrix_B,
            subnetwork,
            common.TopologicalDissimilarityMeasure.BASIC,
        )
        adjusted_regulatory_dissimilarity = (
            subnetwork_statistics.adjusted_regulatory_dissimilarity(
                gene_expression_A,
                gene_expression_B,
                adj_matrix_A,
                adj_matrix_B,
                subnetwork,
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
            "topological_dissimilarity": topological_dissimilarity,
            "adjusted_regulatory_dissimilarity": adjusted_regulatory_dissimilarity,
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
    permutation_test_stats["topological_dissimilarity_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["topological_dissimilarity"]]
        )
    )
    permutation_test_stats["adjusted_regulatory_dissimilarity_norm"] = (
        standardScaler.fit_transform(
            permutation_test_stats[["adjusted_regulatory_dissimilarity"]]
        )
    )

    # Calculate cidrgn
    permutation_test_stats["cidrgn1"] = (
        permutation_test_stats["adjusted_regulatory_dissimilarity_norm"]
        + permutation_test_stats["sam_gs"]
    )
    permutation_test_stats["cidrgn2"] = (
        permutation_test_stats["regulatory_dissimilarity_norm"]
        + permutation_test_stats["topological_dissimilarity_norm"]
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
            "regulatory_dissimilarity",
            "topological_dissimilarity",
            "adjusted_regulatory_dissimilarity",
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
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"
            )

            if simulation_name == "sim_gt-random_subnetworksize-10_numsamples-50":
                continue  # TODO: delete this; already ran this experiment

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


if __name__ == "__main__":
    permutation_tests_simulated_data()
