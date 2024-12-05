import pandas as pd
import numpy as np
import os
import argparse
import common
import r_wrappers

params = {
    "u": 0.0001,  # controls off-diagonal elements of the precision matrix
    "v": 0.9,  # added to the diagonal elements of the precision matrix
    "g": 1,  # for cluster/hub, this is the number of groups; for band this is the bandwidth
}


def generate_subnetwork(n_genes, n_samples, subnetwork_type, folder=None, name=None):
    """
    Generate a common subnetwork. Saves data and plots graph if name/folder provided.
    """
    assert (folder == None) == (
        name == None
    ), "Either both folder and name should be provided or neither should be provided."

    # Generate the common subnetwork
    common_subnetwork = r_wrappers.huge_generator(
        num_samples=n_samples,
        num_variables=n_genes,
        graph_type=subnetwork_type,
        u=params["u"],
        v=params["v"],
        g=params["g"],
    )

    if name is not None:
        # Create a folder if it does not exist
        if not os.path.exists(folder):
            os.makedirs(folder)

        # Save the common subnetwork to an Excel file
        with pd.ExcelWriter(f"{folder}/{name}.xlsx", engine="openpyxl") as writer:
            for sheet_name, data in common_subnetwork.items():
                df = pd.DataFrame(data)
                df.to_excel(writer, sheet_name=sheet_name, index=False, header=False)

        # Save the plot of the undirected graph
        common.plot_undirected_graph(
            common_subnetwork["theta"],
            f"{folder}/{name}.png",
            title="Common Subnetwork",
        )

    return common_subnetwork


# def parse_arguments():
#     parser = argparse.ArgumentParser(description="Generate a simulated graph.")
#     parser.add_argument(
#         "--n_common_subnetwork_genes",
#         type=common.parse_positive_int,
#         required=True,
#         help="Number of genes in each common subnetwork.",
#     )
#     parser.add_argument(
#         "--n_differential_subnetwork_genes",
#         type=common.parse_positive_int,
#         required=True,
#         help="Number of genes in each differential subnetwork.",
#     )
#     parser.add_argument(
#         "--n_common_subnetworks",
#         type=common.parse_positive_int,
#         required=True,
#         help="Number of common subnetworks",
#     )
#     parser.add_argument(
#         "--n_differential_subnetworks",
#         type=common.parse_positive_int,
#         required=True,
#         help="Number of differential subnetworks",
#     )
#     parser.add_argument(
#         "--subnetwork_type",
#         type=common.parse_subnetwork_type,
#         required=True,
#         help="Type of subnetwork  for common subnetworks and differential subnetworks in A",
#     )

#     return parser.parse_args()


def create_simulated_gene_expression_data(
    total_num_genes,
    num_common_subnetwork_genes,
    num_differential_subnetwork_genes,
    num_common_subnetworks,
    num_differential_subnetworks,
    num_samples,
    graph_type,
    simulation_name,
):

    assert (
        total_num_genes
        >= num_common_subnetwork_genes * num_common_subnetworks
        + num_differential_subnetwork_genes * num_differential_subnetworks
    ), "Total number of genes should be greater than the sum of genes in common and differential subnetworks"

    print(f"Generating simulated gene expression data for {simulation_name}")

    simulation_folder = f"data/simulated_gene_expression/{simulation_name}"

    # Set up data frames to hold the gene expression for group A and group B
    column_names = [f"sample_{i}" for i in range(num_samples)]
    column_names.append("subnetwork_id")
    full_gene_expr_A = pd.DataFrame(columns=column_names)
    full_gene_expr_B = pd.DataFrame(columns=column_names)

    # Generate data for 4 common subnetworks, each of 10 genes and 100 samples
    for i in range(1, num_common_subnetworks + 1):
        subnetwork = generate_subnetwork(
            num_common_subnetwork_genes,
            num_samples * 2,
            graph_type,
            folder=simulation_folder,
            name=f"common_subnetwork_{i}",
        )
        subnetwork_gene_expression = pd.DataFrame(subnetwork["data"]).T
        subn_expr_A = subnetwork_gene_expression.iloc[:, :num_samples]
        subn_expr_B = subnetwork_gene_expression.iloc[:, num_samples:]
        subn_expr_A.columns = [f"sample_{j}" for j in range(num_samples)]
        subn_expr_B.columns = [f"sample_{j}" for j in range(num_samples)]
        subn_expr_A["subnetwork_id"] = f"common_subnetwork_{i}"
        subn_expr_B["subnetwork_id"] = f"common_subnetwork_{i}"
        full_gene_expr_A = pd.concat([full_gene_expr_A, subn_expr_A], axis=0)
        full_gene_expr_B = pd.concat([full_gene_expr_B, subn_expr_B], axis=0)

    # Generate data for 6 differential subnetworks, each of 10 genes and 50 samples
    for i in range(1, num_differential_subnetworks + 1):
        # Generate differential subnetwork A (of graph type)
        subnetwork = generate_subnetwork(
            num_differential_subnetwork_genes,
            num_samples,
            graph_type,
            folder=simulation_folder,
            name=f"differential_subnetwork_A_{i}",
        )
        subn_expr_A = pd.DataFrame(subnetwork["data"]).T
        subn_expr_A.columns = [f"sample_{j}" for j in range(num_samples)]
        subn_expr_A["subnetwork_id"] = f"differential_subnetwork_{i}"
        full_gene_expr_A = pd.concat([full_gene_expr_A, subn_expr_A], axis=0)

        # generate differential subnetwork B (of band type)
        subnetwork = generate_subnetwork(
            num_differential_subnetwork_genes,
            num_samples,
            common.SubnetworkType.BAND,
            folder=simulation_folder,
            name=f"differential_subnetwork_B_{i}",
        )
        subn_expr_B = pd.DataFrame(subnetwork["data"]).T
        subn_expr_B.columns = [f"sample_{j}" for j in range(num_samples)]
        subn_expr_B["subnetwork_id"] = f"differential_subnetwork_{i}"
        full_gene_expr_B = pd.concat([full_gene_expr_B, subn_expr_B], axis=0)

    # Sample remaining genes from a normal gaussian distribution
    num_remaining_genes = total_num_genes - (
        num_common_subnetwork_genes * num_common_subnetworks
        + num_differential_subnetwork_genes * num_differential_subnetworks
    )
    remaining_gene_expr_A = pd.DataFrame(
        np.random.normal(size=(num_remaining_genes, num_samples))
    )
    remaining_gene_expr_B = pd.DataFrame(
        np.random.normal(size=(num_remaining_genes, num_samples))
    )
    remaining_gene_expr_A.columns = [f"sample_{j}" for j in range(num_samples)]
    remaining_gene_expr_B.columns = [f"sample_{j}" for j in range(num_samples)]
    remaining_gene_expr_A["subnetwork_id"] = "remaining_genes"
    remaining_gene_expr_B["subnetwork_id"] = "remaining_genes"
    full_gene_expr_A = pd.concat([full_gene_expr_A, remaining_gene_expr_A], axis=0)
    full_gene_expr_B = pd.concat([full_gene_expr_B, remaining_gene_expr_B], axis=0)

    # Organize expression data for entire experiment
    full_gene_expr_A.reset_index(drop=True, inplace=True)
    full_gene_expr_A.insert(
        0, "gene_id", [f"gene_{i}" for i in range(full_gene_expr_A.shape[0])]
    )
    full_gene_expr_B.reset_index(drop=True, inplace=True)
    full_gene_expr_B.insert(
        0, "gene_id", [f"gene_{i}" for i in range(full_gene_expr_B.shape[0])]
    )

    # Create metadata to keep track of gene IDs and subnetwork IDs
    assert full_gene_expr_A["gene_id"].equals(
        full_gene_expr_B["gene_id"]
    ), "Gene IDs do not match between groups A and B"
    assert full_gene_expr_A["subnetwork_id"].equals(
        full_gene_expr_B["subnetwork_id"]
    ), "Subnetwork IDs do not match between groups A and B"
    metadata = pd.DataFrame(
        {
            "gene_id": full_gene_expr_A["gene_id"],
            "subnetwork_id": full_gene_expr_A["subnetwork_id"],
        }
    )

    full_gene_expr_A.drop(columns=["subnetwork_id"], inplace=True)
    full_gene_expr_B.drop(columns=["subnetwork_id"], inplace=True)

    # Save data
    full_gene_expr_A.to_csv(
        f"{simulation_folder}/gene_expression_A.tsv", index=False, sep="\t"
    )
    full_gene_expr_B.to_csv(
        f"{simulation_folder}/gene_expression_B.tsv", index=False, sep="\t"
    )
    metadata.to_csv(f"{simulation_folder}/metadata.tsv", index=False, sep="\t")

    # Plot heatmaps for gene expression data
    full_gene_expr_A.index = full_gene_expr_A["gene_id"]
    full_gene_expr_A.drop(columns=["gene_id"], inplace=True)
    common.plot_gene_expression_heatmap(
        full_gene_expr_A,
        metadata,
        f"{simulation_folder}/gene_expression_A_raw.png",
        title="",
    )
    common.plot_gene_expression_heatmap(
        full_gene_expr_A,
        metadata,
        f"{simulation_folder}/gene_expression_A_standardized.png",
        title="",
        standardized=True,
    )
    full_gene_expr_B.index = full_gene_expr_B["gene_id"]
    full_gene_expr_B.drop(columns=["gene_id"], inplace=True)
    common.plot_gene_expression_heatmap(
        full_gene_expr_B,
        metadata,
        f"{simulation_folder}/gene_expression_B_raw.png",
        title="",
    )
    common.plot_gene_expression_heatmap(
        full_gene_expr_B,
        metadata,
        f"{simulation_folder}/gene_expression_B_standardized.png",
        title="",
        standardized=True,
    )

    # Compute the empirical covariance matrices for the gene expression data
    print("Computing empirical covariance matrices")
    empirical_cov_A = full_gene_expr_A.T.cov()
    empirical_cov_B = full_gene_expr_B.T.cov()
    empirical_cov_A.to_csv(
        f"{simulation_folder}/empirical_cov_A.tsv", index=True, sep="\t"
    )
    empirical_cov_B.to_csv(
        f"{simulation_folder}/empirical_cov_B.tsv", index=True, sep="\t"
    )

    # Plot covariance matrix
    print("Plotting covariance matrices")
    common.plot_covariance_heatmap(
        empirical_cov_A,
        f"{simulation_folder}/empirical_cov_A.png",
        title="Empirical Covariance Matrix A",
    )
    common.plot_covariance_heatmap(
        empirical_cov_B,
        f"{simulation_folder}/empirical_cov_B.png",
        title="Empirical Covariance Matrix B",
    )


def redo_plots(simulated_data_folder):
    metadata = pd.read_csv(f"{simulated_data_folder}/metadata.tsv", sep="\t")

    gene_expression_data_A = pd.read_csv(
        f"{simulated_data_folder}/gene_expression_A.tsv",
        sep="\t",
        index_col=0,
    )
    common.plot_gene_expression_heatmap(
        gene_expression_data_A,
        metadata,
        f"{simulated_data_folder}/gene_expression_heatmap_A.png",
    )

    gene_expression_data_B = pd.read_csv(
        f"{simulated_data_folder}/gene_expression_B.tsv",
        sep="\t",
        index_col=0,
    )
    common.plot_gene_expression_heatmap(
        gene_expression_data_B,
        metadata,
        f"{simulated_data_folder}/gene_expression_heatmap_B.png",
    )

    covariance_matrix_A = pd.read_csv(
        f"{simulated_data_folder}/empirical_cov_A.tsv", sep="\t", index_col=0
    )
    common.plot_covariance_heatmap(
        covariance_matrix_A,
        metadata,
        f"{simulated_data_folder}/empirical_cov_A.png",
    )
    covariance_matrix_B = pd.read_csv(
        f"{simulated_data_folder}/empirical_cov_B.tsv", sep="\t", index_col=0
    )
    common.plot_covariance_heatmap(
        covariance_matrix_B,
        metadata,
        f"{simulated_data_folder}/empirical_cov_B.png",
    )


def create_normal_size_dataset():
    for graph_type in common.SubnetworkType:
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"
            )

            redo_plots(f"data/simulated_gene_expression/{simulation_name}")

            # create_simulated_gene_expression_data(
            #     total_num_genes=1000,
            #     num_common_subnetwork_genes=10,
            #     num_differential_subnetwork_genes=10,
            #     num_common_subnetworks=4,
            #     num_differential_subnetworks=6,
            #     num_samples=50,
            #     graph_type=graph_type,
            #     simulation_name=simulation_name,
            # )


def create_small_dataset(rep=None):
    for graph_type in common.SubnetworkType:
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{4}_numsamples-{20}"
            )
            if rep is not None:
                simulation_name = f"rep{rep}/{simulation_name}"

            create_simulated_gene_expression_data(
                total_num_genes=50,
                num_common_subnetwork_genes=4,
                num_differential_subnetwork_genes=4,
                num_common_subnetworks=4,
                num_differential_subnetworks=0,
                num_samples=20,
                graph_type=graph_type,
                simulation_name=simulation_name,
            )


if __name__ == "__main__":
    # create_small_dataset(rep=5)
    create_normal_size_dataset()
