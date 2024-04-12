import collections
from collections import defaultdict
from typing import Union

import networkx as nx

from utils.common import sort_dict_freq


def sort_node_appearances_in_sub_graph(appearances: list[tuple], neuron_names: list) -> dict[Union[int, str], int]:
    """
    :param appearances: the sub graphs appearances of a given motif
    :param neuron_names: list of neurons names for neural network or an empty list otherwise
    :return: a sorted dict - each key is a node and the value is it frequency
    """
    nodes_count = defaultdict(int)
    for sub_graph in appearances:
        graph = nx.DiGraph()
        graph.add_edges_from(sub_graph)
        for n in list(graph.nodes):
            node = neuron_names[n] if neuron_names else n
            nodes_count[node] += 1

    return sort_dict_freq(nodes_count)


def sort_node_based_on_order_in_sub_graph(appearances: list[tuple],
                                          neuron_names: list,
                                          roles: list[tuple]) -> dict[list]:
    """
    :param appearances: the sub graphs appearances of a given motif
    :param neuron_names: list of neurons names for neural network or an empty list otherwise
    :param roles: list of tuples with the pattern of roles of the motif
    :return: dict, where each key is role, and the value is a sorted list based on appearances of that role
    """
    node_roles = defaultdict(list)
    for sub_graph in appearances:
        neurons_in_sub_graph = {}
        for (src, tar), (role1, role2) in zip(sub_graph, roles):
            neurons_in_sub_graph[role1] = src
            neurons_in_sub_graph[role2] = tar

        for role, n in neurons_in_sub_graph.items():
            node = neuron_names[n] if neuron_names else n
            node_roles[role].append(node)

    freq_node_roles = {}
    for role in node_roles:
        freq_node_roles[role] = sort_dict_freq(dict(collections.Counter(node_roles[role])))

    return freq_node_roles
