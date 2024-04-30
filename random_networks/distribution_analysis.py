import networkx as nx

from networks.loaders.network_loader import NetworkLoader
from utils.config import Config

if __name__ == "__main__":
    config = Config()
    config.set_property("run_args", "plot_properties", "true")
    config.set_property("run_args", "plot_full_graph", "false")
    config.set_property("run_args", "run_sub_graph_search", "false")
    config.set_property("run_args", "run_motif_criteria", "false")

    n = 1000
    e = 4000
    p = e / (n * (n - 1))
    er_graph = nx.erdos_renyi_graph(n=n, p=p, directed=True)
    ba_graph = nx.DiGraph(nx.barabasi_albert_graph(n=n, m=2))  # power law

    loader = NetworkLoader()
    network = loader.load_graph(ba_graph)
