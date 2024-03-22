from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from utils.simple_logger import Logger

sub_graphs_ids = {'feed_forwards': [38, 44, 104, 134, 194, 200]}


def get_sub_id_name(sub_id: int):
    for sub_graph in sub_graphs_ids:
        if sub_id in sub_graphs_ids[sub_graph]:
            return sub_graph
    return None


class SubGraphs(metaclass=ABCMeta):
    def __init__(self, network: DiGraph):
        self.network = network
        self.logger = Logger()

    @abstractmethod
    def search_sub_graphs(self, k: int) -> dict:
        pass
