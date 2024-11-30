import argparse
from enum import Enum
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


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
