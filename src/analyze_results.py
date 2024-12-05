import pandas as pd
import os
import common
from common import SubnetworkType
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

params = {}


def analyze_permutation_test_results(
    permutation_test_stats: pd.DataFrame, p_values: pd.DataFrame
):
    """
    Analyzes the results of the permutation tests

    Parameters:
    - permutation_test_stats: df containing the test statistics for each subnetwork (real and random permutations)
        (indices are subnetwork_id and columns are statistics)
    - p_values: df containing the p-values for each subnetwork (indices are subnetwork_id for subnetworks of interest and
        columns are p-values associated with statistics)

    Returns:
    - data frame with basic analysis results for each stat method
        (indices are methods, columns are results)
    """

    methods = p_values.columns
    analysis_results = pd.DataFrame(index=methods, columns=["TN", "FP", "FN", "TP"])

    p_values["subnetwork_type"] = p_values.index.map(
        lambda x: (
            "common"
            if x.startswith("common")
            else "differential" if x.startswith("differential") else "unknown"
        )
    )

    def count_FP(method):
        return (
            (p_values["subnetwork_type"] == "common")
            & (p_values[method] < 0.05)  # TODO: debug why this isn't working
        ).sum()

    def count_TP(method):
        return (
            (p_values["subnetwork_type"] == "differential") & (p_values[method] < 0.05)
        ).sum()

    def count_FN(method):
        return (
            (p_values["subnetwork_type"] == "differential") & (p_values[method] >= 0.05)
        ).sum()

    def count_TN(method):
        return (
            (p_values["subnetwork_type"] == "common") & (p_values[method] >= 0.05)
        ).sum()

    def auroc(method):
        targets = p_values["subnetwork_type"] == "differential"
        preds = 1 - p_values[method]
        return roc_auc_score(targets, preds)

    analysis_results["TP"] = [count_TP(method) for method in methods]
    analysis_results["FP"] = [count_FP(method) for method in methods]
    analysis_results["FN"] = [count_FN(method) for method in methods]
    analysis_results["TN"] = [count_TN(method) for method in methods]
    analysis_results["AUROC"] = [auroc(method) for method in methods]
    analysis_results["recall"] = analysis_results["TP"] / (
        analysis_results["TP"] + analysis_results["FN"]
    )
    analysis_results["precision"] = analysis_results["TP"] / (
        analysis_results["TP"] + analysis_results["FP"]
    )
    analysis_results["TNR"] = analysis_results["TN"] / (
        analysis_results["TN"] + analysis_results["FP"]
    )
    analysis_results["FPR"] = analysis_results["FP"] / (
        analysis_results["FP"] + analysis_results["TN"]
    )
    analysis_results["Fmeasure"] = (
        2
        * (analysis_results["precision"] * analysis_results["recall"])
        / (analysis_results["precision"] + analysis_results["recall"])
    )
    analysis_results["accuracy"] = (analysis_results["TP"] + analysis_results["TN"]) / (
        analysis_results["TP"]
        + analysis_results["TN"]
        + analysis_results["FP"]
        + analysis_results["FN"]
    )

    return analysis_results


def analyze_results_simulated_data():
    """
    Analyzes the permutation test results for each of the simulated data sets.
    Requires the output from permutation tests to be present
    """
    # Create a multindex df of combined simulated results
    # for comparison with cidrgn table 1 simulated results
    result_stats = ["recall", "precision", "TNR", "Fmeasure", "accuracy", "AUROC"]
    methods = [
        "sam_gs",
        "gsca",
        "cidrgn1",
        "cidrgn2",
        "regulatory_dissimilarity",
        "topological_dissimilarity",
        "adjusted_regulatory_dissimilarity",
    ]
    subnetwork_sizes = [10, 50]
    subnetwork_types = [
        SubnetworkType.RANDOM,
        SubnetworkType.CLUSTER,
        SubnetworkType.SCALE_FREE,
        SubnetworkType.HUB,
    ]
    row_index = pd.MultiIndex.from_product(
        [result_stats, methods], names=["Statistic", "Method"]
    )
    col_index = pd.MultiIndex.from_product(
        [subnetwork_sizes, subnetwork_types],
        names=["Subnetwork size", "Subnetwork Type"],
    )
    combined_results = pd.DataFrame(index=row_index, columns=col_index)

    for subnetwork_type in common.SubnetworkType:
        if subnetwork_type != common.SubnetworkType.BAND:
            simulation_name = f"sim_gt-{subnetwork_type.name.lower()}_subnetworksize-{10}_numsamples-{50}"

            # Read in permutation test results
            permutation_test_stats = pd.read_csv(
                f"data/simulated_data_results/{simulation_name}/permutation_test_stats.tsv",
                sep="\t",
                index_col=0,
            )
            p_values = pd.read_csv(
                f"data/simulated_data_results/{simulation_name}/p_values.tsv",
                sep="\t",
                index_col=0,
            )

            # Perform the analysis
            results = analyze_permutation_test_results(permutation_test_stats, p_values)

            # Save results
            folder = f"data/simulated_data_results/{simulation_name}"
            if not os.path.exists(folder):
                os.makedirs(folder)

            results.to_csv(
                f"{folder}/results.tsv",
                sep="\t",
            )

            # Add to combined results
            for method in methods:
                for result in result_stats:
                    combined_results.loc[(result, method), (10, subnetwork_type)] = (
                        results.loc[method, result]
                    )

    # Save combined results
    combined_results.to_csv(
        f"data/simulated_data_results/combined_results.tsv",
        sep="\t",
    )


if __name__ == "__main__":
    analyze_results_simulated_data()
