import collections

import networkx as nx
from networkx import DiGraph
from subgraphs.sub_graphs_utils import get_sub_graph_mapping_to_motif
from utils.common import get_decimal_from_bin_vec, get_bin_vec_from_decimal
from utils.types import PolarityFrequencies
from collections import defaultdict


def count_network_polarity_ratio(polarity: list[str]):
    polarity_frequencies = collections.Counter(polarity)
    polarity_options = set(polarity_frequencies.keys())
    if polarity_options != {'+', '-'}:
        raise Exception(f'Polarity {str(polarity_options)} not supported!')
    polarity_ratio = polarity_frequencies['+'] / polarity_frequencies['-']
    return polarity_ratio


def get_polarity_frequencies(appearances: list[tuple[tuple]],
                             roles: list[tuple],
                             graph: DiGraph) -> list[PolarityFrequencies]:
    """
    :param appearances: the sub graphs appearances of a given motif
    :param roles: list of tuples with the pattern of roles of the motif
    :param graph: the graph from which the appearances were extracted
    :return: a list of object containing frequencies and polarity list: e.g.: [-, +] for a 2 edge sub graph
    the order in each polarity list is the same as the order of the roles
    """

    edges = len(roles)
    fsl = defaultdict(list)

    for sub_graph in appearances:
        nodes_in_sub_graph = get_sub_graph_mapping_to_motif(sub_graph, roles)
        nodes_in_sub_graph_reverse = {v: k for k, v in nodes_in_sub_graph.items()}

        sub_graph_polarity = {}
        for s, t in sub_graph:
            polarity = graph.get_edge_data(s, t)['polarity']
            role_edge = (nodes_in_sub_graph_reverse[s], nodes_in_sub_graph_reverse[t])
            edge_idx = roles.index(role_edge)  # TODO:  bug in polarity rand network, use m=25 on the pol network
            sub_graph_polarity[edge_idx] = polarity

        polarity_vec = sorted(sub_graph_polarity.values(), key=lambda item: item[0])
        polarity_bin_vec = list(map(lambda x: 1 if x == '+' else 0, polarity_vec))
        polarity_id = get_decimal_from_bin_vec(polarity_bin_vec)
        fsl[polarity_id].append(sub_graph)

    fsl = dict(sorted(fsl.items()))
    polarity_frequencies = []
    for decimal in range((2 ** edges)):
        sub_graphs = fsl.get(decimal, [])
        freq = len(sub_graphs)
        bin_digits = get_bin_vec_from_decimal(decimal=decimal, pad_to=edges)
        polarity_vec = list(map(lambda x: '+' if x == 1 else '-', bin_digits))
        polarity_frequencies.append(PolarityFrequencies(frequency=freq, polarity=polarity_vec, sub_graphs=sub_graphs))

    return polarity_frequencies
