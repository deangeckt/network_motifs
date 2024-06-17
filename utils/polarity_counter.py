import collections
from typing import Union

from isomorphic.isomorphic import get_sub_graph_mapping_to_motif, IsomorphicMotifMatch
from utils.types import PolarityFrequencies
from collections import defaultdict
from itertools import product


def count_network_polarity_ratio(polarity: list[str]) -> collections.Counter:
    polarity_frequencies = collections.Counter(polarity)
    total = sum(polarity_frequencies.values())
    for key in polarity_frequencies:
        polarity_frequencies[key] /= total

    return polarity_frequencies


def get_polarity_frequencies(appearances: list[tuple[tuple]],
                             roles: list[tuple],
                             polarity_options: list[str],
                             motif_id: Union[int, str],
                             iso_matcher: IsomorphicMotifMatch,
                             ) -> list[PolarityFrequencies]:
    """
    :param iso_matcher: an iso matcher object
    :param motif_id: the original motif id (i.e.: not the polarity motif id)
    :param appearances: the sub graphs appearances of a given motif
    :param roles: list of tuples with the pattern of roles of the motif
    :param polarity_options: list of polarity options: can be from [+, -, complex]
    :return: a list of object containing frequencies and polarity list: e.g.: [-, +] for a 2 edge sub graph
    the order in each polarity list is the same as the order of the roles
    """
    edges = len(roles)
    products = list(product(polarity_options, repeat=edges))
    fsl = defaultdict(list)

    for sub_graph in appearances:
        nodes_in_sub_graph = get_sub_graph_mapping_to_motif(sub_graph, roles, [])
        nodes_in_sub_graph_reverse = {v: k for k, v in nodes_in_sub_graph.items()}

        sub_graph_polarity = {}

        for s, t, polarity_attr in sub_graph:
            polarity = polarity_attr['polarity']
            role_edge = (nodes_in_sub_graph_reverse[s], nodes_in_sub_graph_reverse[t])
            edge_idx = roles.index(role_edge)  # TODO: still bug here (try k=4 m=15, [+ -])
            sub_graph_polarity[edge_idx] = polarity

        polarity_vec = tuple([sub_graph_polarity[k] for k in sorted(sub_graph_polarity.keys())])
        polarity_id = products.index(polarity_vec)
        fsl[polarity_id].append(sub_graph)

    polarity_frequencies = []
    for decimal, polarity_vec in enumerate(products):
        sub_graphs = fsl.get(decimal, [])
        freq = len(sub_graphs)
        polarity_frequencies.append(PolarityFrequencies(frequency=freq, polarity=polarity_vec, sub_graphs=sub_graphs))

    return iso_matcher.merge_polarity_isomorphic_frequencies(motif_id, polarity_frequencies)
