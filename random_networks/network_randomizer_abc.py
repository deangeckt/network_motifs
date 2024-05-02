import random
from abc import ABCMeta, abstractmethod

import networkx as nx
from networkx import Graph, DiGraph

from networks.network import Network
from utils.simple_logger import Logger


class NetworkRandomizer(metaclass=ABCMeta):
    def __init__(self, network: Network):
        self.network = network
        self.logger = Logger()

        self.inhibitory_polarity_ratio = 1 / network.polarity_ratio

    @abstractmethod
    def generate(self, amount: int) -> list[DiGraph]:
        pass

    def _log_avg_num_of_generated_edges(self, random_networks: list[DiGraph], amount: int):
        avg_edges = sum([len(rand_network.edges) for rand_network in random_networks]) / amount
        self.logger.info(f'average # edges of all random networks: {round(avg_edges, 3)}')

    @staticmethod
    def _remove_random_edges(graph: Graph, amount: int):
        if amount < 0:
            raise Exception('Cant remove negative amount of edges')

        all_edges = list(graph.edges)
        random.shuffle(all_edges)
        removed_edges = all_edges[:amount]
        graph.remove_edges_from(removed_edges)

    @staticmethod
    def _assign_direction(graph: Graph) -> DiGraph:
        dir_graph = nx.DiGraph()
        for u, v in graph.edges:
            new_edge = (u, v) if random.random() < 0.5 else (v, u)
            dir_graph.add_edge(*new_edge)
        return dir_graph

    def _assign_polarity(self, graph: DiGraph):
        r = self.network.polarity_ratio
        e = len(self.network.graph.edges)

        inhibitory_edges = round(e / (r+1))
        excitatory_edges = e - inhibitory_edges

        all_edges = list(graph.edges)
        random.shuffle(all_edges)

        for ex_idx in range(excitatory_edges):
            edge = all_edges[ex_idx]
            nx.set_edge_attributes(graph, {edge: {"polarity": '+'}})
        for in_idx in range(excitatory_edges, len(all_edges)):
            edge = all_edges[in_idx]
            nx.set_edge_attributes(graph, {edge: {"polarity": '-'}})
