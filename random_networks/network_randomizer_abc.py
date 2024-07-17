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

    @abstractmethod
    def generate(self, amount: int) -> list[DiGraph]:
        pass

    def _log_avg_num_of_generated_edges(self, random_networks: list[DiGraph], amount: int):
        avg_edges = sum([len(rand_network.edges) for rand_network in random_networks]) / amount
        # TODO: check why this is lower than E in pol network... debug with 5 rand networks...
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
        pol_amounts = {}

        curr_sum = 0
        last_key = ''
        for i, key in enumerate(self.network.polarity_ratio):
            last_key = key
            if i == len(self.network.polarity_ratio) - 1:
                break
            pol_amounts[key] = round(self.network.polarity_ratio[key] * len(graph.edges))
            curr_sum += pol_amounts[key]

        pol_amounts[last_key] = len(graph.edges) - curr_sum

        all_edges = list(graph.edges)
        random.shuffle(all_edges)

        start_idx = 0
        for pol in pol_amounts:
            amount = pol_amounts[pol]
            edges = all_edges[start_idx: start_idx + amount]
            edges_dict = {e: {"polarity": pol} for e in edges}
            nx.set_edge_attributes(graph, edges_dict)
            start_idx += amount
