from abc import ABCMeta, abstractmethod

from networkx import DiGraph

from networks.network import Network
from utils.simple_logger import Logger


class NetworkRandomizer(metaclass=ABCMeta):
    def __init__(self, network: Network):
        self.network = network
        self.logger = Logger()

    @abstractmethod
    def generate(self, amount: int) -> list[DiGraph]:
        pass
