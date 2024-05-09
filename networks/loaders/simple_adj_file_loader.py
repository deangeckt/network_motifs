from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from networks.network import Network


class SimpleAdjFileLoader(NetworkLoaderStrategy):
    def __init__(self, args):
        super().__init__(args)

    def load(self, *args) -> Network:
        """
        simple txt format: (v1, v2, w) per line
        whereas v1 -> v2, and w is ignored
        """
        file_path = args[0]
        with open(file_path, "r") as f:
            for line in f.readlines():
                v1, v2, _ = tuple(line.strip().split())
                self.graph.add_edge(int(v1), int(v2))

        return self._copy_network_params()
