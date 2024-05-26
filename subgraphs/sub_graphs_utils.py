from collections import defaultdict

import networkx as nx
import numpy as np
from networkx import DiGraph

from tqdm import tqdm

from utils.common import get_decimal_from_bin_vec, get_bin_vec_from_decimal
from utils.types import MotifName, Motif
from networkx.algorithms import isomorphism


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


two_sub_graphs_ids = {
    MotifName.mutual_regulation: [6],
}

three_sub_graphs_ids = {
    MotifName.feed_forward: [38, 44, 104, 134, 194, 200],
    MotifName.cascade: [12, 34, 66, 96, 132, 136],
    MotifName.fan_out: [6, 40, 192],
    MotifName.fan_in: [36, 72, 130]
}

four_sub_graphs_ids = {
    MotifName.bi_fan: [204, 2448, 2570, 13056, 20560, 24582],
    MotifName.bi_parallel: [904, 2136, 2182, 4544, 6672, 8716, 10498, 12356, 16458, 16532, 20994, 24848],
    MotifName.sim_3: [14, 208, 2816, 28672]
}

sub_graphs_ids_per_k = {
    2: two_sub_graphs_ids,
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


def get_adjacency_matrix(graph: DiGraph) -> np.ndarray:
    nodes = list(graph.nodes)
    nodes.sort()
    return nx.adjacency_matrix(graph, nodelist=nodes).todense()


def get_id(graph: DiGraph) -> int:
    adj_mat = get_adjacency_matrix(graph)
    vec = adj_mat.flatten()
    return get_decimal_from_bin_vec(list(vec))


def get_sub_graph_from_id(decimal: int, k: int) -> DiGraph:
    bin_digits = get_bin_vec_from_decimal(decimal=decimal, pad_to=k ** 2)
    bin_digits.reverse()
    adj_mat = np.array(bin_digits).reshape(k, k)
    return nx.DiGraph(adj_mat)


def generate_isomorphic_k_sub_graphs(k: int, allow_self_loops=False) -> tuple[dict, dict]:
    """
    :param allow_self_loops: edges from a node to itself
    :param k: motif / sub graph size
    :return: isomorphic_mapping: a dict where each key is a motif / sub graph id and
    the value is the smallest motif id of the same isomorphic set of sub graphs.
    :return isomorphic_graphs: a dict where each key is the smallest motif id and the value
    is the list of all the motif ids that are isomorphic to it.

    not scalable for k >= 5.
    In case we are NOT looking for Anti Motifs - this isn't mandatory for mfinder.
    """
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


def get_number_of_disjoint_group_nodes(sub_graphs: list[tuple[tuple]]) -> int:
    """
    :param sub_graphs: all the sub graphs (list of edges) participating in a candidate motif sub graph
    :return: the number of completely disjoint groups of nodes
    """
    graph = nx.DiGraph()
    for sub_graph in sub_graphs:
        graph.add_edges_from(sub_graph)

    return nx.number_weakly_connected_components(graph)


def get_role_pattern(adj_mat: np.ndarray) -> list[tuple]:
    """
    :param adj_mat: the adjacency matrix of the motif
    :return: list of tuples with roles of the motif, in the format: (a,b), (b,c)
    """
    ascii_start = 97
    roles: list[tuple] = []
    for src, arr in enumerate(adj_mat):
        for tar, val in enumerate(arr):
            if val != 1:
                continue
            roles.append((chr(src + ascii_start), chr(tar + ascii_start)))
    return roles


def get_sub_graph_mapping_to_motif(sub_graph: tuple[tuple], motif_roles: list[tuple]) -> dict:
    """
    :param sub_graph: tuple of tuples - the edges in the subgraph
    :param motif_roles: list of tuples with roles of a motif, in the format: (a,b), (b,c)
    :return: mapping of each role per node: e.g: {'a' : node_id, 'b' : node_id}
    """
    graph = nx.DiGraph()
    graph.add_edges_from(sub_graph)

    motif_graph = nx.DiGraph()
    motif_graph.add_edges_from(motif_roles)

    matcher = isomorphism.GraphMatcher(motif_graph, graph)
    if not (matcher.is_isomorphic()):
        raise Exception('The sub graph is not isomorphic to the motif')
    return dict(matcher.mapping)


def create_base_motif(sub_id: int, k: int) -> Motif:
    name = get_sub_id_name(sub_id=sub_id, k=k)
    sub_graph = get_sub_graph_from_id(decimal=sub_id, k=k)
    adj_mat = nx.adjacency_matrix(sub_graph).todense()
    role_pattern = get_role_pattern(adj_mat)
    return Motif(name=name, id=sub_id, adj_mat=adj_mat, role_pattern=role_pattern)


def create_sim_motif(sim_id: str, adj_mat: np.ndarray):
    roles = get_role_pattern(adj_mat)
    # we don't add the adj_mat, so it won't explode the table in the logs
    return Motif(name=MotifName.na, id=sim_id, adj_mat=np.array([]), role_pattern=roles)
