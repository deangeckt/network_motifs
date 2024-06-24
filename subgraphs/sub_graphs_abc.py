from abc import ABCMeta, abstractmethod
from collections import defaultdict

from networkx import DiGraph
import networkx as nx

from utils.sub_graphs import get_id, get_sub_graph_from_id
from utils.simple_logger import Logger
from utils.types import SubGraphSearchResult


class SubGraphsABC(metaclass=ABCMeta):
    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        self.network = network
        s, t = list(network.edges)[0]
        self.use_polarity = 'polarity' in network[s][t] and network[s][t]['polarity'] is not None

        self.isomorphic_mapping = isomorphic_mapping
        self.logger = Logger()

        self.k = -1  # motif size
        self.allow_self_loops = False

        self.fsl = defaultdict(int)
        self.fsl_fully_mapped = defaultdict(list)

        self.inc_canonical_label_foo = self.__inc_count_w_canonical_label_using_iso_mapping if isomorphic_mapping else \
            self.__inc_count_w_canonical_label_self_iso

    @abstractmethod
    def search_sub_graphs(self, k: int, allow_self_loops: bool) -> SubGraphSearchResult:
        """
        :param k: motif size
        :param allow_self_loops: allow or not. effects the self isomorphic version of canonical labeling
        :return: SubGraphSearchResult
        """
        pass

    def _inc_count_w_canonical_label(self, sub_graph: DiGraph):
        self.inc_canonical_label_foo(sub_graph)

    def __inc_count_w_canonical_label_self_iso(self, sub_graph: DiGraph):
        if list(nx.selfloop_edges(sub_graph)) and not self.allow_self_loops:
            return

        sub_id = get_id(sub_graph)
        if len(self.fsl) == 0:
            self.fsl[sub_id] = 1
            self.__append_to_fully_mapped_fsl(sub_id, sub_graph)
            return

        if sub_id in self.fsl:
            self.fsl[sub_id] += 1
            self.__append_to_fully_mapped_fsl(sub_id, sub_graph)
            return

        found = False
        for id_ in self.fsl:
            other_graph = get_sub_graph_from_id(id_, k=self.k)
            if nx.is_isomorphic(sub_graph, other_graph):
                self.fsl[id_] += 1
                self.__append_to_fully_mapped_fsl(id_, sub_graph)
                found = True
                break

        if not found:
            self.fsl[sub_id] = 1
            self.__append_to_fully_mapped_fsl(sub_id, sub_graph)

    def __append_to_fully_mapped_fsl(self, sub_id_isomorphic_representative: int, sub_graph: DiGraph):
        if not self.use_polarity:
            self.fsl_fully_mapped[sub_id_isomorphic_representative].append(tuple(list(sub_graph.edges)))
        else:
            polarities = nx.get_edge_attributes(sub_graph, 'polarity')
            pol_edges = tuple([(*e, {'polarity': polarities[e]}) for e in list(sub_graph.edges)])
            self.fsl_fully_mapped[sub_id_isomorphic_representative].append(pol_edges)

    def __inc_count_w_canonical_label_using_iso_mapping(self, sub_graph: DiGraph):
        sub_id = get_id(sub_graph)
        if sub_id not in self.isomorphic_mapping:
            return

        sub_id_isomorphic_representative = self.isomorphic_mapping[sub_id]
        self.fsl[sub_id_isomorphic_representative] += 1
        self.__append_to_fully_mapped_fsl(sub_id_isomorphic_representative, sub_graph)
