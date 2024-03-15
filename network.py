import numpy as np


class Network:
    def __init__(self):
        self.nodes = None
        self.adj_matrix = None

        # debug: in case nodes names aren't arranged serially
        self.mapped_nodes_reverse = None
        self.mapped_nodes = None

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
