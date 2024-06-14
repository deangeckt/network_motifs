import collections
from collections import defaultdict
from typing import Union

import networkx as nx
from networkx.algorithms import isomorphism
from tqdm import tqdm

from utils.sub_graphs import get_sub_graph_from_id
from utils.types import Motif, PolarityFrequencies


def generate_isomorphic_k_sub_graphs(k: int, allow_self_loops=False) -> tuple[dict, dict]:
    """
    :param allow_self_loops: edges from a node to itself
    :param k: motif / sub graph size
    :return: isomorphic_mapping: a dict where each key is a motif / sub graph id and
    the value is the smallest motif id of the same isomorphic set of sub graphs.
    :return isomorphic_graphs: a dict where each key is the smallest motif id and the value
    is the list of all the motif ids that are isomorphic to it.

    not scalable for k >= 5.
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


def get_fsl_ids_iso_mapping(src_fsl_ids: list[int], tar_fsl_ids: list[int], k: int) -> dict[int, int]:
    """
    Used when use_isomorphic_mapping flag is off, i.e.: larger K's.

    :param src_fsl_ids: a list of sub graphs id's of an FSL list
    :param tar_fsl_ids: a list of sub graphs id's of an FSL list
    :param k: motif / sub graph size
    :return: a dict mapping the ids in the source list to the target based on iso matching.
    in case a missing in the target list, the value will be -1
    """
    ids_iso_mappings = {}
    for src in src_fsl_ids:
        if src in tar_fsl_ids:
            ids_iso_mappings[src] = src
            continue

        ids_iso_mappings[src] = -1
        src_graph = get_sub_graph_from_id(decimal=src, k=k)
        for tar in tar_fsl_ids:
            tar_graph = get_sub_graph_from_id(decimal=tar, k=k)
            if nx.is_isomorphic(src_graph, tar_graph):
                ids_iso_mappings[src] = tar
                break
    return ids_iso_mappings


def _get_polarity_sub_graph(role_pattern: list[tuple], polarity: list[str]) -> nx.DiGraph:
    graph = nx.DiGraph(role_pattern)
    for role, pol in zip(role_pattern, polarity):
        s, t = role
        graph[s][t]['polarity'] = pol
    return graph


def get_sub_graph_mapping_to_motif(sub_graph: tuple[tuple],
                                   motif_roles: list[tuple],
                                   polarity: list[str],
                                   ) -> dict:
    """
    :param sub_graph: tuple of tuples - the edges in the subgraph
    :param motif_roles: list of tuples with roles of a motif, in the format: (a,b), (b,c)
    :param polarity: for polarity motifs, the polarity list of string representation
    :return: mapping of each role per node: e.g: {'a' : node_id, 'b' : node_id}
    """
    graph = nx.DiGraph(sub_graph)

    if polarity:
        motif_graph = _get_polarity_sub_graph(motif_roles, polarity)
        # TODO: support fast sim matching?
        matcher = isomorphism.DiGraphMatcher(motif_graph, graph,
                                             edge_match=lambda e1, e2: e1['polarity'] == e2['polarity'])
    else:
        motif_graph = nx.DiGraph(motif_roles)
        matcher = isomorphism.DiGraphMatcher(motif_graph, graph)

    if not (matcher.is_isomorphic()):
        raise Exception('The sub graph is not isomorphic to the motif')
    return dict(matcher.mapping)


def _is_polarity_sim_isomorphic(g1, g2):
    """
    check if two graphs of a SIM polarity motif are isomorphic.
    returns True if number of symbols, e.g.: #"+" == #"-". otherwise False.
    """
    g1_pol_counter = collections.Counter(nx.get_edge_attributes(g1, 'polarity').values())
    g2_pol_counter = collections.Counter(nx.get_edge_attributes(g2, 'polarity').values())
    return g1_pol_counter == g2_pol_counter


def merge_polarity_isomorphic_sub_graphs(
        motif_id: Union[int, str],
        roles: list[tuple],
        polarity_frequencies: list[PolarityFrequencies]) -> list[PolarityFrequencies]:
    isomorphic = []
    del_idx_list = []

    for curr_idx, pol_obj in enumerate(polarity_frequencies):
        sub_graph = _get_polarity_sub_graph(roles, pol_obj.polarity)

        found = False
        found_idx = -1
        for found_idx, reps_graph in isomorphic:
            if isinstance(motif_id, str):
                if _is_polarity_sim_isomorphic(reps_graph, sub_graph):
                    found = True
                    break
            else:
                if nx.is_isomorphic(reps_graph, sub_graph, edge_match=lambda e1, e2: e1['polarity'] == e2['polarity']):
                    found = True
                    break

        if found:
            merge_to_pol_obj = polarity_frequencies[found_idx]
            merge_to_pol_obj.frequency += pol_obj.frequency
            merge_to_pol_obj.sub_graphs.extend(pol_obj.sub_graphs)
            del_idx_list.append(curr_idx)

        if not found:
            isomorphic.append((curr_idx, sub_graph))

    polarity_frequencies = [polarity_frequencies[i] for i in range(len(polarity_frequencies)) if i not in del_idx_list]
    return polarity_frequencies
