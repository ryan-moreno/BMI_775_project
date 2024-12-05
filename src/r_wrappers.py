import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import numpy as np
import pandas as pd

# Activate pandas2ri conversion
pandas2ri.activate()

required_packages = ["huge"]

# Check packages are installed, if not install them
for package in required_packages:
    if not rpackages.isinstalled(package):
        ro.r(f'install.packages("{package}")')

# Import packages
packages = {}
for package in required_packages:
    packages[package] = rpackages.importr(package)

# GSCA must be loaded from the source RData file
ro.r("load('references/GSCA.RData')")
# gsca = ro.r["singleDC"]

# SAM_GS must be loaded from the source R file
ro.r("source('references/sam_gs.R')")
sam_gs_r = ro.r["SAMGS"]


def huge_generator(
    num_samples,
    num_variables,
    graph_type,
    v=None,
    u=None,
    g=1,
):
    """
    Wrapper for the R function huge.generator in Python, with added parameters v, u, and vis.

    Parameters:
    - num_samples: Number of samples (rows)
    - num_variables: Number of variables (columns)
    - graph_type: The graph structure (enum)
    - v: The off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with u. The default value is 0.3.
    - u: A positive number being added to the diagonal elements of the precision matrix, to control the magnitude of partial correlations. The default value is 0.1.
    - g: The number of groups for cluster/hub, or the bandwidth for band. The default value is 1.

    Returns:
    - A dictionary containing the following:
        - data: The n by d matrix for the generated data
        - sigma: The covariance matrix for the generated data
        - omega: The precision matrix for the generated data
        - sigmahat: The empirical covariance matrix for the generated data
        - theta: The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data
    """

    # Get expected string name for graph_type
    graph_type = graph_type.name.lower().replace("_", "-")
    assert graph_type in [
        "random",
        "hub",
        "cluster",
        "band",
        "scale-free",
    ], "Invalid graph_type. Must be one of 'random', 'hub', 'cluster', 'band', or 'scale-free'."

    # Convert required parameters num_samples, num_variables
    r_num_samples = ro.IntVector([num_samples])
    r_num_variables = ro.IntVector([num_variables])

    # Convert optional parameters v, u, vis
    if v is not None:
        r_v = ro.FloatVector([v])
    else:
        r_v = None

    if u is not None:
        r_u = ro.FloatVector([u])
    else:
        r_u = None

    # Call huge.generator in R
    result = packages["huge"].huge_generator(
        r_num_samples, r_num_variables, graph_type, v=r_v, u=r_u, g=1, prob=1
    )

    # Extract elements from the result (S3 object "sim")
    data = np.array(result[0])
    sigma = np.array(result[1])
    omega = np.array(result[2])
    sigmahat = np.array(result[3])
    theta = ro.r["as.matrix"](result[4])

    # Return the results as a dictionary
    output = {
        "data": data,
        "sigma": sigma,
        "omega": omega,
        "sigmahat": sigmahat,
        "theta": theta,
    }

    return output


# def single_study_GSCA(gene_expr_A, gene_expr_B, subnetwork):
#     """
#     Wrapper for the R function single.study.gsca in Python.
#     TODO: this doesn't work because the package doesn't include C symbol dEuc2sampleperm

#     Parameters:
#     - gene_expr_A: Gene expression data for group A
#     - gene_expr_B: Gene expression data for group B
#     - subnetwork: list of gene ids in the subnetwork

#     Returns:
#     - GSCA dispersion index
#     """

#     # Format gene expression data
#     data = gene_expr_A.merge(gene_expr_B)
#     group = [gene_expr_A.shape[1], gene_expr_B.shape[1]]
#     gene_set_definitions = [subnetwork]

#     # Convert data to R matrix
#     r_data = ro.r["as.matrix"](pandas2ri.py2rpy(data))

#     # Convert lists to R vectors
#     r_group = ro.IntVector(group)
#     r_gene_set_definitions = ro.ListVector(
#         {
#             f"gene_set_{i+1}": ro.StrVector(inner_list)
#             for i, inner_list in enumerate(gene_set_definitions)
#         }
#     )

#     # Call singleDC in R
#     result = gsca(r_data, r_group, r_gene_set_definitions, 1)

#     # Extract elements from the result
#     dispersion_index = result[0][0]
#     p_value = result[1][0]
#     permutation_DI_matrix = result[2]

#     # Return dispersion index
#     return dispersion_index


def load_gsca_example_data():
    """
    Load example data for GSCA from the RData file.
    Using expected dispersion indices printed in
    https://www.biostat.wisc.edu/~kendzior/GSCA/GSCA_vignette.pdf

    Returns:
    - A dictionary containing the following:
        - gene_expr_A: Gene expression data for group A
        - gene_expr_B: Gene expression data for group B
        - gene_sets: list of lists of gene ids
        - expected_stats: list of expected dispersion indices
    """

    # Load example data from RData file
    ro.r("load('references/GSCA.RData')")

    # Extract gene expression data from R environment
    r_data = ro.r("LungCancer3$data$Michigan")
    colnames = ro.r("colnames(LungCancer3$data$Michigan)")
    rownames = ro.r("rownames(LungCancer3$data$Michigan)")
    data = pd.DataFrame(r_data, columns=colnames, index=rownames)
    gene_expr_A = data.filter(like="Normal")
    gene_expr_B = data.filter(like="Tumor")

    # Extract first 3 gene sets
    r_gene_sets = ro.r("LungCancer3$info$GSdef")
    gene_sets = [[str(gene) for gene in gene_set] for gene_set in r_gene_sets]
    gene_sets = gene_sets[:5]

    # Expected first 5 dispersion indices
    expected_stats = [0.431, 0.434, 0.404, 0.427, 0.422]

    # Return the data as a dictionary
    output = {
        "gene_expr_A": gene_expr_A,
        "gene_expr_B": gene_expr_B,
        "gene_sets": gene_sets[2:5],  # These are the ones that are under 500 genes
        "expected_stats": expected_stats[2:5],
    }

    return output


def sam_gs(gene_expr_A, gene_expr_B, gene_set):
    """
    Calls the author's R implementation of SAM_GS

    Parameters:
    - gene_expr_A: Gene expression data for group A
    - gene_expr_B: Gene expression data for group B
    - gene_sets: list of lists of gene ids

    Returns:
    - list of SAM_GS scores
    """
    gene_expr_A.columns = [f"A_{str(col)}" for col in gene_expr_A.columns]
    gene_expr_B.columns = [f"B_{str(col)}" for col in gene_expr_B.columns]
    full_expression_data = pd.concat([gene_expr_A, gene_expr_B], axis=1)
    group_labels = ["1"] * gene_expr_A.shape[1] + ["2"] * gene_expr_B.shape[1]

    # Convert data to R
    expr_r = ro.pandas2ri.py2rpy(full_expression_data)
    group_labels_r = ro.StrVector(group_labels)
    gene_sets_r = ro.ListVector(
        {f"gene_set_{i+1}": ro.StrVector(genes) for i, genes in enumerate([gene_set])}
    )

    # Call SAM_GS in R
    result = sam_gs_r(
        gene_sets_r,
        expr_r,
        group_labels_r,
    )

    # Extract t_like_stat per gene from the result
    t_like_stat_per_gene = result[2][0]

    # Sum over all genes to get the total SAM_GS score
    full_expression_data["t_like_stat_gene"] = t_like_stat_per_gene
    full_expression_data["t_like_squared"] = (
        full_expression_data["t_like_stat_gene"] ** 2
    )
    sam_gs_score = full_expression_data[full_expression_data.index.isin(gene_set)][
        "t_like_squared"
    ].sum()

    return sam_gs_score
