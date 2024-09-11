from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.types import NetworkLoaderArgs


class SimpleAdjFileLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        super().__init__(args)

    def load(self, *args):
        """
        simple txt format: (v1, v2, w) per line
        where v1 -> v2, and w is ignored
        """
        file_path = args[0]
        with open(file_path, "r") as f:
            for line in f.readlines():
                v1, v2, _ = tuple(line.strip().split())
                self.network.graph.add_edge(int(v1), int(v2))
