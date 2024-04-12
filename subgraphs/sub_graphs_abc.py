from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from utils.simple_logger import Logger
from utils.types import SubGraphSearchResult


class SubGraphsABC(metaclass=ABCMeta):
    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        self.network = network
        self.isomorphic_mapping = isomorphic_mapping
        self.logger = Logger()

    @abstractmethod
    def search_sub_graphs(self, k: int) -> SubGraphSearchResult:
        """
        :param k: motif size
        :return: SubGraphSearchResult
        """
        pass
