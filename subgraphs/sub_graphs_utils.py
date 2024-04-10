from collections import defaultdict

import networkx as nx
import numpy as np
from networkx import DiGraph

from tqdm import tqdm

from utils.config import Config
from utils.types import MotifName, Motif


class HashedGraph:
    def __init__(self, graph: tuple):
        self.graph = graph

    def __eq__(self, other):
        return sorted(self.graph) == sorted(other.graph)

    def __hash__(self):
        return hash(tuple(sorted(self.graph)))

    def __str__(self):
        return str(self.graph)


def graph_to_hashed_graph(g: DiGraph) -> HashedGraph:
    graph_tuple = tuple(list(g.edges))
    return HashedGraph(graph_tuple)


three_sub_graphs_ids = {
    MotifName.feed_forward: [38, 44, 104, 134, 194, 200],
    MotifName.cascade: [12, 34, 66, 96, 132, 136],
    MotifName.fan_out: [6, 40, 192],
    MotifName.fan_in: []  # TODO: add 36?
}

four_sub_graphs_ids = {
    MotifName.bi_fan: [204, 2448, 2570, 13056, 20560, 24582]
}

sub_graphs_ids_per_k = {
    3: three_sub_graphs_ids,
    4: four_sub_graphs_ids
}


def get_sub_id_name(sub_id: int, k: int) -> MotifName:
    if k not in sub_graphs_ids_per_k:
        return MotifName.na
    sub_graphs_id = sub_graphs_ids_per_k[k]
    for sub_graph_name in sub_graphs_id:
        if sub_id in sub_graphs_id[sub_graph_name]:
            return sub_graph_name
    return MotifName.na


def get_id(graph: DiGraph) -> int:
    nodes = list(graph.nodes)
    nodes.sort()
    adj_mat = nx.adjacency_matrix(graph, nodelist=nodes).todense()
    vec = adj_mat.flatten()
    decimal = 0
    for i, bit in enumerate(vec):
        decimal += bit * (2 ** i)
    return decimal


def get_sub_graph_from_id(decimal: int, k: int) -> DiGraph:
    bin_digits = [int(d) for d in str(bin(decimal))[2:]]
    pad_amount = k ** 2 - len(bin_digits)
    padding = pad_amount * [0]
    bin_digits = padding + bin_digits
    bin_digits.reverse()
    adj_mat = np.array(bin_digits).reshape(k, k)
    return nx.DiGraph(adj_mat)


def generate_isomorphic_k_sub_graphs(k: int) -> tuple[dict, dict]:
    """
    :param k: motif / sub graph size
    :return: isomorphic_mapping: a dict where each key is a motif / sub graph id and
    the value is the smallest motif id of the same isomorphic set of sub graphs.
    :return isomorphic_graphs: a dict where each key is the smallest motif id and the value
    is the list of all the motif ids that are isomorphic to it.

    not scalable for k >= 5.
    In case we are NOT looking for Anti Motifs - this isn't mandatory for mfinder.
    """
    config = Config()
    allow_self_loops = config.get_boolean_property('run_args', 'allow_self_loops')

    isomorphic = defaultdict(dict)

    possible_options = (2 ** (k ** 2))
    for sub_id in tqdm(range(possible_options)):
        sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
        un_dir_ub_graph = nx.Graph(sub_graph)

        if list(nx.selfloop_edges(sub_graph)) and not allow_self_loops:
            continue

        # remove the not connected cases, e.g.: no edges at all-sub graph
        if not nx.is_connected(un_dir_ub_graph):
            continue

        found = False
        for reps_graph_id in isomorphic:
            reps_graph = isomorphic[reps_graph_id]['reps_graph']
            if nx.is_isomorphic(reps_graph, sub_graph):
                isomorphic[reps_graph_id]['list'].append(sub_id)
                found = True
                break

        if not found:
            isomorphic[sub_id]['reps_graph'] = sub_graph
            isomorphic[sub_id]['list'] = [sub_id]

    isomorphic_mapping = {}
    isomorphic_graphs = {}
    for id_ in isomorphic:
        isomorphic_graphs[id_] = isomorphic[id_]['list']
        for sub_id in isomorphic[id_]['list']:
            isomorphic_mapping[sub_id] = id_

    return isomorphic_mapping, isomorphic_graphs


def get_number_of_disjoint_group_nodes(sub_graphs: list[list[tuple]]) -> int:
    """
    :param sub_graphs: all the sub graphs (list of edges) participating in a candidate motif sub graph
    :return: the number of completely disjoint groups of nodes
    """
    graph = nx.DiGraph()
    for sub_graph in sub_graphs:
        graph.add_edges_from(sub_graph)

    return nx.number_weakly_connected_components(graph)


def create_base_motif(sub_id: int, k: int) -> Motif:
    name = get_sub_id_name(sub_id=sub_id, k=k)
    sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
    adj_mat = str(nx.adjacency_matrix(sub_graph).todense())
    return Motif(name=name, id=sub_id, adj_mat=adj_mat)
