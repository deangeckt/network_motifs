from networkx import DiGraph

from subgraphs.sub_graphs_abc import SubGraphsABC
import networkx as nx

from utils.sub_graphs import create_base_motif
from utils.types import SubGraphSearchResult
import netsci.metrics.motifs as nsm
import itertools


# https://github.com/gialdetti/netsci
class NetsciWrapper(SubGraphsABC):
    """
    This class wrap Netsci package and uses the louzoun algorithm which support accessing the actual sub graphs.
    This class support k=3 only without self loops!

    netsci package:
    * Gal, E., Perin, R., Markram, H., London, M., and Segev, I. (2019). Neuron Geometry Underlies a Universal Local Architecture in Neuronal Networks. BioRxiv 656058.

    louzoun algorithm:
    * Itzhack, R., Mogilevski, Y., and Louzoun, Y. (2007). An optimal algorithm for counting network motifs.
      Phys. A Stat. Mech. Its Appl. 381, 482-490.
    """

    def __init__(self, network: DiGraph, isomorphic_mapping: dict):
        super().__init__(network, isomorphic_mapping)
        # self.network.remove_edges_from(nx.selfloop_edges(network))
        # TODO: a bug when input graph has a self loop...

        self.motif_keys = [12, 36, 6, 38, 14, 74, 98, 78, 102, 46, 108, 110, 238]
        self.unique_roles = {'a', 'b', 'c'}

    def _get_actual_edges(self, roles: list[tuple], nodes: list[int]):
        # Generate all possible mappings of roles to nodes
        role_permutations = []
        for perm in itertools.permutations(nodes, len(self.unique_roles)):
            role_map = {role: node for role, node in zip(sorted(self.unique_roles), perm)}
            role_permutations.append(role_map)

        # For each mapping, check if it produces edges that exist in the graph
        for role_map in role_permutations:
            candidate_edges = []
            valid_mapping = True

            for src_role, dst_role in roles:
                src_node = role_map[src_role]
                dst_node = role_map[dst_role]

                # Check if this edge exists in the graph
                if not self.network.has_edge(src_node, dst_node):
                    valid_mapping = False
                    break

                candidate_edges.append((src_node, dst_node))

            if valid_mapping:
                return candidate_edges

        # If no valid mapping is found
        return []

    def search_sub_graphs(self, k: int, allow_self_loops: bool) -> SubGraphSearchResult:
        if k != 3:
            raise Exception('Netsci Wrapper support k=3 only')
        if allow_self_loops:
            raise Exception('Netsci Wrapper does not support self loops')

        sorted_nodes = list(self.network.nodes)
        sorted_nodes.sort()
        A = nx.adjacency_matrix(self.network, nodelist=sorted_nodes).todense()

        n_reals, participating_nodes = nsm.motifs(A, algorithm='louzoun', participation=True)
        n_reals = n_reals[3:]
        participating_nodes = participating_nodes[3:]

        fsl_fully_mapped = {}
        for i, sub_graphs in enumerate(participating_nodes):
            motif_id = self.motif_keys[i]
            motif = create_base_motif(sub_id=motif_id, k=3)
            role_pattern = motif.role_pattern

            fsl_fully_mapped[motif_id] = []
            for nodes in sub_graphs:
                graph_nodes = [sorted_nodes[n] for n in nodes]
                sub_graph_edges = self._get_actual_edges(role_pattern, graph_nodes)
                fsl_fully_mapped[motif_id].append(sub_graph_edges)

        fsl = {self.motif_keys[i]: amount for (i, amount) in enumerate(n_reals)}
        return SubGraphSearchResult(fsl=fsl, fsl_fully_mapped=fsl_fully_mapped)
