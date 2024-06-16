import collections
from collections import defaultdict
from typing import Union

import networkx as nx

from utils.common import sort_dict_freq
from isomorphic.isomorphic import get_sub_graph_mapping_to_motif
from utils.types import Motif


def sort_node_appearances_in_sub_graph(appearances: list[tuple[tuple]],
                                       neuron_names: list) -> dict[Union[int, str], int]:
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


def sort_node_roles_in_sub_graph(appearances: list[tuple[tuple]],
                                 neuron_names: list,
                                 motif: Motif
                                 ) -> dict[str, dict]:
    """
    :param appearances: the sub graphs appearances of a given motif
    :param neuron_names: list of neurons names for neural network or an empty list otherwise
    :param motif: motif object with roles: list of tuples with the pattern of roles of the motif
    :return: dict, where each key is role, and the value is a sorted dict based on appearances of that role
    """
    node_roles = defaultdict(list)
    for sub_graph in appearances:
        nodes_in_sub_graph = get_sub_graph_mapping_to_motif(sub_graph, motif.role_pattern, motif.polarity)
        for role, n in nodes_in_sub_graph.items():
            node = neuron_names[n] if neuron_names else n
            node_roles[role].append(node)

    freq_node_roles = {}
    for role in node_roles:
        freq_node_roles[str(role)] = sort_dict_freq(dict(collections.Counter(node_roles[role])))

    return freq_node_roles

