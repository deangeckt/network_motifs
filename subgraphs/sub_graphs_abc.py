from abc import ABCMeta, abstractmethod
from collections import defaultdict

from networkx import DiGraph
import networkx as nx

from subgraphs.sub_graphs_utils import get_id
from utils.simple_logger import Logger
from utils.types import SubGraphSearchResult


class SubGraphsABC(metaclass=ABCMeta):
    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        self.network = network
        s, t = list(network.edges)[0]
        self.use_polarity = network[s][t]['polarity'] is not None

        self.isomorphic_mapping = isomorphic_mapping
        self.logger = Logger()
        self.k = -1  # motif size

        self.fsl = defaultdict(int)
        self.fsl_fully_mapped = defaultdict(list)

    @abstractmethod
    def search_sub_graphs(self, k: int) -> SubGraphSearchResult:
        """
        :param k: motif size
        :return: SubGraphSearchResult
        """
        pass

    def _inc_count_w_canonical_label(self, sub_graph: DiGraph):
        sub_id = get_id(sub_graph)
        if sub_id not in self.isomorphic_mapping:
            return

        sub_id_isomorphic_representative = self.isomorphic_mapping[sub_id]
        self.fsl[sub_id_isomorphic_representative] += 1

        if not self.use_polarity:
            self.fsl_fully_mapped[sub_id_isomorphic_representative].append(tuple(list(sub_graph.edges)))
        else:
            polarities = nx.get_edge_attributes(sub_graph, 'polarity')
            pol_edges = tuple([(*e, {'polarity': polarities[e]}) for e in list(sub_graph.edges)])
            self.fsl_fully_mapped[sub_id_isomorphic_representative].append(pol_edges)