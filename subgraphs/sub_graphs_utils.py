from collections import defaultdict
from typing import Optional, Tuple, Dict, Any

import networkx as nx
import numpy as np
from networkx import DiGraph
from enum import Enum

from tqdm import tqdm

from utils.config import Config


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


class SubGraphAlgoName(str, Enum):
    specific = 'specific'
    mfinder_induced = 'mfinder_i'
    mfinder_none_induced = 'mfinder_ni'


class MotifName(str, Enum):
    self_loops = 'self_loops'
    mutual_regulation = 'mutual_regulation'
    fan_outs = 'fan_outs'
    cascades = 'cascades'
    feed_forwards = 'feed_forwards'
    bi_fan = 'bi_fan'


three_sub_graphs_ids = {
    MotifName.feed_forwards: [38, 44, 104, 134, 194, 200],
    MotifName.cascades: [12, 34, 66, 96, 132, 136],
    MotifName.fan_outs: [6, 40, 192]
}

four_sub_graphs_ids = {
    MotifName.bi_fan: [204, 2448, 2570, 13056, 20560, 24582]
}

sub_graphs_ids_per_k = {
    3: three_sub_graphs_ids,
    4: four_sub_graphs_ids
}


def get_sub_id_name(sub_id: int, k: int) -> Optional[MotifName]:
    if k not in sub_graphs_ids_per_k:
        return None
    sub_graphs_id = sub_graphs_ids_per_k[k]
    for sub_graph_name in sub_graphs_id:
        if sub_id in sub_graphs_id[sub_graph_name]:
            return sub_graph_name
    return None


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

    # TODO: uniQ still not aligned with mfinder software

    def sub_graph_to_nodes_set(sub_graph: list[tuple]) -> set:
        nodes = set([v1 for v1, v2 in sub_graph])
        nodes.update(set([v2 for v1, v2 in sub_graph]))
        return nodes

    uniq = 0
    all_sets: list[set] = [sub_graph_to_nodes_set(sub_graph) for sub_graph in sub_graphs]

    for sub_set in all_sets:
        all_others_as_a_set = set()
        for set_ in all_sets:
            if set_ == sub_set:
                continue
            all_others_as_a_set.update(set_)

        for node in sub_set:
            if node not in all_others_as_a_set:
                uniq += 1
                print(f'Uniq node: {node} in the sub_set: {sub_set}')
                break

    return uniq
