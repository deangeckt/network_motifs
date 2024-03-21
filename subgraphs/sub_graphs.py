from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from utils.simple_logger import Logger


class SubGraphs(metaclass=ABCMeta):
    def __init__(self, network: DiGraph):
        self.network = network
        self.logger = Logger()

    @abstractmethod
    def search_sub_graphs(self, k: int) -> dict:
        pass
