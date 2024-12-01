# import pandas as pd
import numpy as np
import os
from sklearn.linear_model import Lasso
import common
import pandas as pd

params = {
    "regularization_param": 0.5,  # regularization parameter for LASSO based on simulated data
}


def estimate_regularization_parameter():
    """
    Estimate the regularization parameter for LASSO.
    """

    trials_list = []

    # Find optimal regularization parameter for each simulated dataset
    for graph_type in common.SubnetworkType:
        if graph_type != common.SubnetworkType.BAND:
            simulation_name = (
                f"sim_gt-{graph_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"
            )

            # Read in metadata
            metadata = pd.read_csv(
                f"data/simulated_gene_expression/{simulation_name}/metadata.tsv",
                sep="\t",
            )
            subnetwork_ids = metadata["subnetwork_id"].unique()
            subnetwork_ids = [
                subnetwork_id
                for subnetwork_id in subnetwork_ids
                if subnetwork_id != "remaining_genes"
            ]

            # GRN on A and B phenotypes separately
            for phenotype in ["A", "B"]:
                gene_expr = pd.read_csv(
                    f"data/simulated_gene_expression/{simulation_name}/gene_expression_{phenotype}.tsv",
                    sep="\t",
                    index_col="gene_id",
                )

                # Get the number of edges in the true regulatory network
                true_number_edges = 0
                for subnetwork_id in subnetwork_ids:
                    if "differential" in subnetwork_id:
                        split_subnetwork_id = subnetwork_id.rsplit("_", 1)
                        subnetwork_string = f"{split_subnetwork_id[0]}_{phenotype}_{split_subnetwork_id[1]}"
                    else:
                        subnetwork_string = subnetwork_id

                    true_subnetwork_adj_matrix = pd.read_excel(
                        f"data/simulated_gene_expression/{simulation_name}/{subnetwork_string}.xlsx",
                        sheet_name="theta",
                    )

                    true_number_edges += np.count_nonzero(true_subnetwork_adj_matrix)

                # Perform LASSO
                for regularization_param in np.linspace(0.1, 1, 10):
                    predicted_adj_matrix = infer_grn_lasso(
                        gene_expr, regularization_param
                    )
                    predicted_num_edges = np.count_nonzero(predicted_adj_matrix)
                    trials_list.append(
                        {
                            "simulation_name": simulation_name,
                            "phenotype": phenotype,
                            "regularization_param": regularization_param,
                            "predicted_num_edges": predicted_num_edges,
                            "true_num_edges": true_number_edges,
                        },
                    )

    # Save results
    trials_df = pd.DataFrame(trials_list)
    path = "data/simulated_inferred_networks/lasso_regularization_param_trials.tsv"
    trials_df.to_csv(path, sep="\t")


def infer_grn_simulated_data():
    for graph_type in common.SubnetworkType:
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

            # Perform LASSO
            adj_matrix_A = infer_grn_lasso(gene_expr_A, params["regularization_param"])
            adj_matrix_B = infer_grn_lasso(gene_expr_B, params["regularization_param"])

            # Save results
            folder = f"data/simulated_inferred_networks/{simulation_name}"
            if not os.path.exists(folder):
                os.makedirs(folder)

            adj_matrix_A.to_csv(
                f"{folder}/adj_matrix_A.tsv",
                sep="\t",
            )
            adj_matrix_B.to_csv(
                f"{folder}/adj_matrix_B.tsv",
                sep="\t",
            )


def infer_grn_lasso(gene_expression, regularization_param):
    """
    Infer the gene regulatory network using LASSO.
    TODO: This does not currently separate target and regulator genes

    Parameters:
    - gene_expression: pandas DataFrame containing gene expression data (rows = genes, columns = samples)
    - regularization_param: float, regularization parameter for LASSO

    Returns:
    - adj_matrix: pandas DataFrame containing the inferred gene regulatory network coefficients
        (value at (i, j) is the coefficient for the effect of gene j on gene i)
    """

    assert regularization_param >= 0, "Regularization parameter must be non-negative"

    regulatory_network = pd.DataFrame(
        np.zeros((gene_expression.shape[0], gene_expression.shape[0])),
        index=gene_expression.index,
        columns=gene_expression.index,
    )
    for target_gene_id in gene_expression.index:
        target_gene_expr = gene_expression.loc[target_gene_id, :]
        target_gene_idx = gene_expression.index.get_loc(target_gene_id)
        predictors = np.delete(gene_expression.values, target_gene_idx, axis=0).T

        # Fit LASSO model
        lasso = Lasso(alpha=regularization_param)
        lasso.fit(predictors, target_gene_expr)

        # Get coefficients from LASSO
        coefficients = np.insert(lasso.coef_, target_gene_idx, 0)
        regulatory_network.loc[target_gene_id, :] = coefficients

    num_edges = np.count_nonzero(regulatory_network)
    print(
        f"Number of edges in predicted regulatory network: {num_edges}/{gene_expression.shape[0]**2}"
    )

    return regulatory_network


if __name__ == "__main__":
    infer_grn_simulated_data()
    # estimate_regularization_parameter()
