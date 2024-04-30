import networkx as nx

from networks.loaders.network_loader import NetworkLoader
from utils.config import Config

if __name__ == "__main__":
    loader = NetworkLoader()

    config = Config()
    config.set_property("run_args", "plot_properties", "true")
    config.set_property("run_args", "plot_full_graph", "false")
    config.set_property("run_args", "run_sub_graph_search", "false")
    config.set_property("run_args", "run_motif_criteria", "false")

    n = 1000
    e = 4000
    p = e / (n * (n - 1))
    er_graph = nx.erdos_renyi_graph(n=n, p=p, directed=True)  # random
    ws_graph = nx.DiGraph(nx.watts_strogatz_graph(n=n, k=10, p=0.5))  # small world, but not scale free
    ba_graph = nx.DiGraph(nx.barabasi_albert_graph(n=n, m=2))  # scale free, i.e.: power law

    print('erdos_renyi')
    loader.load_graph(er_graph)
    print('watts_strogatz')
    loader.load_graph(ws_graph)
    print('barabasi_albert')
    loader.load_graph(ba_graph)

