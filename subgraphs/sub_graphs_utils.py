from typing import Optional

import networkx as nx
import numpy as np
from networkx import DiGraph
from enum import Enum


class SubGraphAlgo(str, Enum):
    specific = 'specific'
    mfinder = 'mfinder'


class MotifName(str, Enum):
    self_loops = 'self_loops'
    mutual_regulation = 'mutual_regulation'
    fan_outs = 'fan_outs'
    cascades = 'cascades'
    feed_forwards = 'feed_forwards'


three_sub_graphs_ids = {MotifName.feed_forwards: [38, 44, 104, 134, 194, 200]}


# def rename_sub_graphs(decimal_dict: dict, k: int) -> dict:
#     """
#     rename sub id's to motif names if exist
#     """
#     res = {}
#     for decimal_id in decimal_dict:
#         sub_name = get_sub_id_name(decimal_id, k)
#         sub_key = sub_name if sub_name else decimal_id
#         res[sub_key] = decimal_dict[decimal_id]
#     return res


def get_sub_id_name(sub_id: int, k: int) -> Optional[MotifName]:
    if k != 3:
        return None
    for sub_graph_name in three_sub_graphs_ids:
        if sub_id in three_sub_graphs_ids[sub_graph_name]:
            return sub_graph_name
    return None


def get_id(sub_graph: tuple) -> int:
    graph = nx.DiGraph(list(sub_graph))
    adj_mat = nx.adjacency_matrix(graph).todense()
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
