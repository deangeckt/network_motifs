from typing import Optional

import networkx as nx
import numpy as np
from networkx import DiGraph
from enum import Enum


class UniqueSubGraph:
    def __init__(self, sub_graph: tuple):
        self.sub_graph = sub_graph

    def __eq__(self, other):
        return sorted(self.sub_graph) == sorted(other.sub_graph)

    def __hash__(self):
        return hash(tuple(sorted(self.sub_graph)))

    def __str__(self):
        return str(self.sub_graph)


class SubGraphAlgoName(str, Enum):
    specific = 'specific'
    mfinder = 'mfinder'


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
    MotifName.bi_fan: [204]  # TODO: add the rest
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


def get_id(sub_graph: tuple) -> int:
    graph = nx.DiGraph(list(sub_graph))
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
