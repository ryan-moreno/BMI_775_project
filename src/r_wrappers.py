import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
import numpy as np

# Check if the Huge package is installed, and if not, install it
if not rpackages.isinstalled("huge"):
    ro.r('install.packages("huge")')

# Import Huge R package
huge = rpackages.importr("huge")


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
    result = huge.huge_generator(
        r_num_samples, r_num_variables, graph_type, v=r_v, u=r_u, g=1
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
