from network import Network
from utils.simple_logger import Logger


class SubGraphs:
    def __init__(self, network: Network, use_logger=True):
        self.network = network
        self.logger = Logger()
        self.logger.toggle(use_logger)

    