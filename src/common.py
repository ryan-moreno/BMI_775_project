import argparse
from enum import Enum
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns


class SubnetworkType(Enum):
    RANDOM = 1
    CLUSTER = 2
    SCALE_FREE = 3
    HUB = 4
    BAND = 5


class StatisticMethod(Enum):
    GSCA = 1
    SAM_GS = 2
    CIDRGN1 = 3
    CIDRGN2 = 4


class TopologicalDissimilarityMeasure(Enum):
    BASIC = 1


def parse_positive_int(value):
    try:
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
        return ivalue
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")


def parse_subnetwork_type(value):
    try:
        # convert input to enum
        return SubnetworkType(value.lower())
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid SubnetworkType '{value}', must be one of {', '.join([e.name for e in SubnetworkType])}"
        )


def plot_undirected_graph(adj_matrix, path, title=""):
    """
    Plots an undirected graph using networkx from a given adjacency matrix.

    Parameters:
    - adj_matrix: A 2D numpy array representing the adjacency matrix.
    - title: Title for the plot (default no title).
    """
    graph = nx.from_numpy_array(np.array(adj_matrix))

    plt.figure(figsize=(8, 6))
    nx.draw(
        graph,
        with_labels=False,
        node_size=500,
        node_color="skyblue",
        edge_color="gray",
    )

    plt.title(title)
    plt.savefig(path)
    plt.close()


def plot_directed_graph(adj_matrix, path, node_groups, title=""):
    """
    Plots a directed graph using networkx from a given adjacency matrix.
    Colors nodes and edges by common and differential subnetworks.

    Parameters:
    - adj_matrix: A 2D numpy array representing the adjacency matrix.
    - node_groups: A pandas DataFrame with columns 'gene_id', 'subnetwork_id'.
    - title: Title for the plot (default no title).
    """
    graph = nx.from_pandas_adjacency(adj_matrix, create_using=nx.DiGraph)
    nx.set_node_attributes(
        graph,
        {row.gene_id: row.subnetwork_id for _, row in node_groups.iterrows()},
        "subnetwork_id",
    )
    nx.set_node_attributes(
        graph,
        {row.gene_id: 1 for _, row in node_groups.iterrows()},
        "level",  # Want to plot the grey nodes underneath the colored nodes
    )

    if node_groups is not None:
        # Create a color map; gray for background genes
        unique_group_ids = node_groups["subnetwork_id"].unique()
        color_map = {subnetwork_id: "gray" for subnetwork_id in unique_group_ids}

        # Warm colors for common subnetworks
        cmap = plt.get_cmap("autumn")
        ids_for_shared_networks = [
            id for id in unique_group_ids if id.startswith("common")
        ]
        color_map_shared = {
            subnetwork_id: cmap(i / len(ids_for_shared_networks))
            for i, subnetwork_id in enumerate(ids_for_shared_networks)
        }
        color_map.update(color_map_shared)

        # Cool colors for differential subnetworks
        cmap = plt.get_cmap("winter")
        ids_for_differential_networks = [
            id for id in unique_group_ids if id.startswith("differential")
        ]
        color_map_differential = {
            subnetwork_id: cmap(i / len(ids_for_differential_networks))
            for i, subnetwork_id in enumerate(ids_for_differential_networks)
        }
        color_map.update(color_map_differential)

        nx.set_node_attributes(
            graph,
            {
                row.gene_id: color_map[row.subnetwork_id]
                for _, row in node_groups.iterrows()
            },
            "color",
        )

        # Edges will be grey if they connect two non-subnetwork nodes or two subnetwork nodes
        # of different subnetworks; otherwise they will be colored like the subnetwork they touch
        edge_colors = {}
        for e in graph.edges():
            u, v = e
            u_color = graph.nodes[u]["color"]
            v_color = graph.nodes[v]["color"]
            if u_color == v_color:
                edge_colors[e] = "grey"
            else:
                if u_color == "grey":
                    edge_colors[e] = v_color
                elif v_color == "grey":
                    edge_colors[e] = u_color
                else:
                    edge_colors[e] = "grey"

        nx.set_edge_attributes(graph, edge_colors, "color")

    # Order nodes so that the grey nodes are plotted first
    node_order = sorted(
        graph.nodes(),
        key=lambda x: 0 if graph.nodes[x]["color"] == "gray" else 1,
    )
    node_colors = [graph.nodes[node]["color"] for node in node_order]

    # Order edges so that the grey edges are plotted first
    edge_order = sorted(
        graph.edges(),
        key=lambda x: 0 if graph.edges[x]["color"] == "grey" else 1,
    )
    edge_colors = [graph.edges[edge]["color"] for edge in edge_order]

    plt.figure(figsize=(8, 6))
    nx.draw(
        graph,
        with_labels=False,
        node_size=50,
        node_color=node_colors,
        edge_color=edge_colors,
        nodelist=node_order,
        edgelist=edge_order,
    )

    plt.title(title)
    plt.savefig(path)
    plt.close()


def plot_covariance_heatmap(covariance_matrix, metadata, path, title=""):
    """
    Plots a heatmap of a covariance matrix.

    Parameters:
    - covariance_matrix: A 2D numpy array representing the covariance matrix.
    - metadata: A pandas DataFrame with columns 'gene_id', 'subnetwork_id'.
    """

    # Cluster genes by their subnetwork ID
    metadata.index = metadata["gene_id"]
    covariance_matrix = covariance_matrix[metadata.sort_values("subnetwork_id").index]

    # Color genes by their subnetwork ID
    unique_subnetwork_ids = metadata["subnetwork_id"].unique()
    cmap = plt.get_cmap("tab20")
    color_map = {
        subnetwork_id: cmap(i / len(unique_subnetwork_ids))
        for i, subnetwork_id in enumerate(unique_subnetwork_ids)
    }
    gene_colors = covariance_matrix.columns.map(
        lambda x: color_map[metadata.loc[x, "subnetwork_id"]]
    )

    # Plot the heatmap
    plt.figure(figsize=(8, 6))
    plot = sns.clustermap(
        covariance_matrix,
        cmap="coolwarm",
        annot=False,
        xticklabels=False,
        yticklabels=False,
        col_colors=gene_colors,
        row_cluster=False,
        col_cluster=False,
    )
    plot.ax_heatmap.set_ylabel("Genes")
    plot.ax_heatmap.set_xlabel("Genes")

    plt.title(title)
    plt.savefig(path)
    plt.close()


def plot_gene_expression_heatmap(
    gene_expression_data,
    metadata,
    path,
    title="",
    standardized=False,
):
    """
    Plots a heatmap of gene expression data.

    Parameters:
    - gene_expression_data: A 2D numpy array representing gene expression data.
    - metadata: A pandas DataFrame with columns 'gene_id', 'subnetwork_id'.
    - path: Path to save file
    - title: Title for the plot (default no title).
    - standardized: Whether to standardize the data before plotting (default False).
    """

    gene_expression_data = gene_expression_data.T

    if standardized:
        # Change data to be std dev from mean so it's easier to compare covariance
        gene_expression_data = gene_expression_data.apply(
            lambda x: (x - x.mean()) / x.std()
        )

    # Cluster genes by their subnetwork ID
    metadata.index = metadata["gene_id"]
    gene_expression_data = gene_expression_data[
        metadata.sort_values("subnetwork_id").index
    ]

    # Color genes by their subnetwork ID
    unique_subnetwork_ids = metadata["subnetwork_id"].unique()
    cmap = plt.get_cmap("tab20")
    color_map = {
        subnetwork_id: cmap(i / len(unique_subnetwork_ids))
        for i, subnetwork_id in enumerate(unique_subnetwork_ids)
    }
    gene_colors = gene_expression_data.columns.map(
        lambda x: color_map[metadata.loc[x, "subnetwork_id"]]
    )

    plot = sns.clustermap(
        gene_expression_data,
        col_cluster=False,  # For debugging
        row_cluster=False,
        cmap="coolwarm",
        figsize=(12, 12),
        annot=False,
        col_colors=gene_colors,
        xticklabels=False,
        yticklabels=False,
    )
    plot.ax_heatmap.set_ylabel("Samples")
    plot.ax_heatmap.set_xlabel("Genes")

    # Legend for subnetwork IDs
    handles = [
        mpatches.Patch(color=color, label=subnetwork)
        for subnetwork, color in color_map.items()
    ]
    plt.legend(handles=handles, title="Subnetwork ID", bbox_to_anchor=(3.7, 0))

    plt.title(title)
    plt.savefig(path)
    plt.close()
