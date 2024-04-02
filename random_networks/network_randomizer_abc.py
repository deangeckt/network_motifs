from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from utils.simple_logger import Logger


class NetworkRandomizer(metaclass=ABCMeta):
    def __init__(self, network: DiGraph):
        self.network = network
        self.logger = Logger()

    @abstractmethod
    def generate(self, amount: int) -> list[DiGraph]:
        pass
