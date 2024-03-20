import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from simple_logger import Logger


class Network:
    def __init__(self):
        self.logger = Logger()

        self.nodes: list[list[int]] = []
        self.adj_matrix: np.ndarray = np.empty(0)

        # debug: in case nodes names aren't arranged serially
        self.mapped_nodes_reverse: dict = {}
        self.mapped_nodes: dict = {}

    def load_graph(self, txt: list[str]):
        nodes_set = set()
        for line in txt:
            v1, v2, w = tuple(line.strip().split())
            nodes_set.add(int(v1))
            nodes_set.add(int(v2))

        sorted_nodes = list(nodes_set)
        sorted_nodes.sort()
        self.mapped_nodes = {n: idx for idx, n in enumerate(sorted_nodes)}
        self.mapped_nodes_reverse = {idx: n for idx, n in enumerate(sorted_nodes)}

        N = len(nodes_set)
        self.adj_matrix = np.zeros((N, N))
        self.nodes = [[] for _ in range(N)]

        for line in txt:
            v1, v2, w = tuple(line.strip().split())
            v1_idx = self.mapped_nodes[int(v1)]
            v2_idx = self.mapped_nodes[int(v2)]

            self.adj_matrix[v1_idx, v2_idx] = 1
            self.nodes[v1_idx].append(v2_idx)

    def log_properties(self):
        self.logger.info(f'Network properties:')
        self.logger.info(f'  - Nodes: {len(self.nodes)}')
        self.logger.info(f'  - Edges: {np.count_nonzero(self.adj_matrix)}')
        self.logger.info('')

    def get_edges(self) -> list:
        """
        return a list of tuples with the graph edges
        """
        edges = []
        for i, node in enumerate(self.nodes):
            for neighbor in node:
                edges.append((self.mapped_nodes_reverse[i], self.mapped_nodes_reverse[neighbor]))
        return edges

    def plot(self):
        G = nx.DiGraph()
        G.add_edges_from(self.get_edges())
        nx.draw_networkx(G, node_size=600, node_color='lightgreen')
        plt.show()
