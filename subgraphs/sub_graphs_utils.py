import networkx as nx
import numpy as np
from networkx import DiGraph

three_sub_graphs_ids = {'feed_forwards': [38, 44, 104, 134, 194, 200]}


def get_sub_id_name(sub_id: int, k: int):
    if k != 3:
        return None
    for sub_graph in three_sub_graphs_ids:
        if sub_id in three_sub_graphs_ids[sub_graph]:
            return sub_graph
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
