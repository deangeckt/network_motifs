import networkx as nx
from scipy import sparse

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.types import NetworkLoaderArgs


class ScipySparseLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        Loads a custom intersected network
        """
        super().__init__(args)

    def load(self, *args):
        file_path = args[0]
        matrix = sparse.load_npz(file_path)

        self.network.graph = nx.from_scipy_sparse_array(matrix, create_using=nx.DiGraph)
