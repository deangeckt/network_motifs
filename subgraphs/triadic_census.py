from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from utils.types import SubGraphSearchResult


class TriadicCensus(SubGraphsABC):
    """
    triadic_census
    Vladimir Batagelj and Andrej Mrvar, A subquadratic triad census algorithm for large sparse networks with small maximum degree, University of Ljubljana
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        self.triadic_key_to_motif_id = {
            '021D': 6,
            '021U': 36,
            '021C': 12,
            '111D': 74,
            '111U': 14,
            '030T': 38,
            '030C': 98,
            '201': 78,
            '120D': 108,
            '120U': 46,
            '120C': 102,
            '210': 110,
            '300': 238
        }

    def search_sub_graphs(self, k: int) -> SubGraphSearchResult:
        if k != 3:
            raise Exception('Triadic Census support k=3 only')

        triadic_census = nx.triadic_census(self.network)
        del triadic_census['003']
        del triadic_census['012']
        del triadic_census['102']

        fsl = dict((self.triadic_key_to_motif_id[key], value) for (key, value) in triadic_census.items())
        return SubGraphSearchResult(fsl=fsl, fsl_fully_mapped={})
