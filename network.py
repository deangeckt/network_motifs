import networkx as nx
import matplotlib.pyplot as plt
from utils.simple_logger import Logger


class Network:
    def __init__(self):
        self.logger = Logger()
        self.graph = nx.DiGraph()

    def load_from_txt(self, txt: list[str]):
        for line in txt:
            v1, v2, w = tuple(line.strip().split())
            self.graph.add_edge(int(v1), int(v2))

    def log_properties(self):
        self.logger.info(f'Network properties:')
        self.logger.info(f'  - Nodes: {len(self.graph)}')
        self.logger.info(f'  - Edges: {len(self.graph.edges)}')
        self.logger.info('')

    def plot(self):
        nx.draw_networkx(self.graph, with_labels=True, node_size=600, node_color='lightgreen')
        plt.show()
