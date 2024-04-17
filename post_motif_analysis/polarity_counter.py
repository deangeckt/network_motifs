import networkx as nx
from networkx import DiGraph
from subgraphs.sub_graphs_utils import get_sub_graph_mapping_to_motif
from utils.common import get_decimal_from_bin_vec, get_bin_vec_from_decimal
from utils.types import PolarityFrequencies


def get_polarity_frequencies(appearances: list[tuple[tuple]],
                             roles: list[tuple],
                             graph: DiGraph) -> list[PolarityFrequencies]:
    """
    :param appearances: the sub graphs appearances of a given motif
    :param roles: list of tuples with the pattern of roles of the motif
    :param graph: the real network graph
    :return: a list of object containing frequencies and polarity list: e.g.: [-, +] for a 2 edge sub graph
    the order in each polarity list is the same as the order of the roles
    """

    edges = len(roles)
    freq_list = [0] * (2 ** edges)

    for sub_graph in appearances:
        nodes_in_sub_graph = get_sub_graph_mapping_to_motif(sub_graph, roles)
        nodes_in_sub_graph_reverse = {v: k for k, v in nodes_in_sub_graph.items()}

        sub_graph_polarity = {}
        for s, t in sub_graph:
            polarity = graph.get_edge_data(s, t)['polarity']
            role_edge = (nodes_in_sub_graph_reverse[s], nodes_in_sub_graph_reverse[t])
            edge_idx = roles.index(role_edge)
            sub_graph_polarity[edge_idx] = polarity

        polarity_vec = list(sub_graph_polarity.values())
        polarity_bin_vec = list(map(lambda x: 1 if x == '+' else 0, polarity_vec))
        polarity_id = get_decimal_from_bin_vec(polarity_bin_vec)
        freq_list[polarity_id] += 1

    polarity_frequencies = []
    for decimal, freq in enumerate(freq_list):
        bin_digits = get_bin_vec_from_decimal(decimal=decimal, pad_to=edges)
        polarity_vec = list(map(lambda x: '+' if x == 1 else '-', bin_digits))
        polarity_frequencies.append(PolarityFrequencies(frequency=freq, polarity=polarity_vec))

    return polarity_frequencies


def get_all_sub_graph_polarities(sub_graphs: list[tuple[tuple]],
                                 graph: nx.DiGraph) -> list:
    """
    attaching the polarity to each edge in the sub graph
    """
    enriched_sub_graphs = []

    for sub_graph in sub_graphs:
        enriched_sub_graph = tuple()
        for edge in sub_graph:
            s, t = edge
            polarity = graph.get_edge_data(s, t)['polarity']
            polarity = '+' if polarity == 1 else '-'
            edge += (polarity,)
            enriched_sub_graph += (edge,)
        enriched_sub_graphs.append(enriched_sub_graph)

    return enriched_sub_graphs
