from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from utils.simple_logger import Logger


class SubGraphsABC(metaclass=ABCMeta):
    def __init__(self, network: DiGraph):
        self.network = network
        self.logger = Logger()

    @abstractmethod
    def search_sub_graphs(self, k: int) -> dict:
        """
        :param k: motif size
        :return: a dict where each key is a sub graph id and value is the frequency of that sub graph
        """
        pass

    @abstractmethod
    def get_sub_graphs_fully_mapped(self) -> dict:
        """
        :return: a dict where each key is a sub graph id and value a list with all sub graphs (list of tuple of edges)
        """
        pass

