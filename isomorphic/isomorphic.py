import collections
import pickle
from collections import defaultdict
from itertools import product
from typing import Union, TypedDict

import networkx as nx
from networkx.algorithms import isomorphism
from tqdm import tqdm

from utils.sub_graphs import get_sub_graph_from_id, create_base_motif
from utils.types import PolarityFrequencies


class IsomorphicMappingBinFile(TypedDict):
    isomorphic_mapping: dict
    isomorphic_graphs: dict


class PolarityIsoMappingBinFile(TypedDict):
    mapping: dict


def get_fsl_ids_iso_mapping(src_fsl_ids: list[int], tar_fsl_ids: list[int], k: int) -> dict[int, int]:
    """
    Used when use_isomorphic_mapping flag is off, i.e.: larger K's.
    #TODO: does this support polarity large k? do we even want this anymore...?
    :param src_fsl_ids: a list of sub graphs id's of an FSL list
    :param tar_fsl_ids: a list of sub graphs id's of an FSL list
    :param k: motif / sub graph size
    :return: a dict mapping the ids in the source list to the target based on iso matching.
    in case a missing in the target list, the value will be -1
    """
    ids_iso_mappings = {}
    for src in src_fsl_ids:
        if src in tar_fsl_ids:
            ids_iso_mappings[src] = src
            continue

        ids_iso_mappings[src] = -1
        src_graph = get_sub_graph_from_id(decimal=src, k=k)
        for tar in tar_fsl_ids:
            tar_graph = get_sub_graph_from_id(decimal=tar, k=k)
            if nx.is_isomorphic(src_graph, tar_graph):
                ids_iso_mappings[src] = tar
                break
    return ids_iso_mappings


def _get_polarity_sub_graph(role_pattern: list[tuple], polarity: list[str]) -> nx.DiGraph:
    graph = nx.DiGraph(role_pattern)
    for role, pol in zip(role_pattern, polarity):
        s, t = role
        graph[s][t]['polarity'] = pol
    return graph


def get_sub_graph_mapping_to_motif(sub_graph: tuple[tuple],
                                   motif_roles: list[tuple],
                                   polarity: list[str],
                                   ) -> dict:
    """
    :param sub_graph: tuple of tuples - the edges in the subgraph
    :param motif_roles: list of tuples with roles of a motif, in the format: (a,b), (b,c)
    :param polarity: for polarity motifs, the polarity list of string representation
    :return: mapping of each role per node: e.g: {'a' : node_id, 'b' : node_id}
    """
    graph = nx.DiGraph(sub_graph)

    if polarity:
        motif_graph = _get_polarity_sub_graph(motif_roles, polarity)
        matcher = isomorphism.DiGraphMatcher(motif_graph, graph,
                                             edge_match=lambda e1, e2: e1['polarity'] == e2['polarity'])
    else:
        motif_graph = nx.DiGraph(motif_roles)
        matcher = isomorphism.DiGraphMatcher(motif_graph, graph)

    if not (matcher.is_isomorphic()):
        raise Exception('The sub graph is not isomorphic to the motif')
    return dict(matcher.mapping)


def _is_polarity_sim_isomorphic(g1, g2):
    """
    check if two graphs of a SIM polarity motif are isomorphic. FASTER than the standard checking
    returns True if number of symbols are equal. e.g.: #"+" == #"-". otherwise False.
    """
    g1_pol_counter = collections.Counter(nx.get_edge_attributes(g1, 'polarity').values())
    g2_pol_counter = collections.Counter(nx.get_edge_attributes(g2, 'polarity').values())
    return g1_pol_counter == g2_pol_counter


class IsomorphicMotifMatch:
    def __init__(self, k: int, polarity_options: list[str], allow_self_loops=False):
        self.base_path = 'isomorphic/mapping'
        self.k = k
        self.allow_self_loops = allow_self_loops
        self.polarity_options = polarity_options

        # Load iso mapping
        file_path = self.__get_isomorphic_k_file_name()
        bin_file: IsomorphicMappingBinFile = self.__import_iso_mapping(file_path)

        self.isomorphic_mapping = bin_file['isomorphic_mapping']
        self.isomorphic_graphs = bin_file['isomorphic_graphs']

        # Load polarity iso mapping
        if not polarity_options:
            return

        file_path = self.__get_polarity_isomorphic_k_file_name()
        bin_file: PolarityIsoMappingBinFile = self.__import_pol_iso_mapping(file_path)
        self.polarity_mapping = bin_file['mapping']

    def generate_isomorphic_k_sub_graphs(self):
        """
        :export: isomorphic_mapping: a dict where each key is a motif / sub graph id and
        the value is the smallest motif id of the same isomorphic set of sub graphs.
        :export isomorphic_graphs: a dict where each key is the smallest motif id and the value
        is the list of all the motif ids that are isomorphic to it.

        not scalable for k >= 5.
        """
        isomorphic = defaultdict(dict)

        possible_options = (2 ** (self.k ** 2))
        for sub_id in tqdm(range(possible_options)):
            sub_graph = get_sub_graph_from_id(decimal=sub_id, k=self.k)
            un_dir_ub_graph = nx.Graph(sub_graph)

            if list(nx.selfloop_edges(sub_graph)) and not self.allow_self_loops:
                continue

            # remove the not connected cases, e.g.: no edges at all-sub graph
            if not nx.is_connected(un_dir_ub_graph):
                continue

            found = False
            for reps_graph_id in isomorphic:
                reps_graph = isomorphic[reps_graph_id]['reps_graph']
                if nx.is_isomorphic(reps_graph, sub_graph):
                    isomorphic[reps_graph_id]['list'].append(sub_id)
                    found = True
                    break

            if not found:
                isomorphic[sub_id]['reps_graph'] = sub_graph
                isomorphic[sub_id]['list'] = [sub_id]

        isomorphic_mapping = {}
        isomorphic_graphs = {}
        for id_ in isomorphic:
            isomorphic_graphs[id_] = isomorphic[id_]['list']
            for sub_id in isomorphic[id_]['list']:
                isomorphic_mapping[sub_id] = id_

        file_path = self.__get_isomorphic_k_file_name()
        self.__export_iso_mapping(file_path, IsomorphicMappingBinFile(isomorphic_mapping=isomorphic_mapping,
                                                                      isomorphic_graphs=isomorphic_graphs))

    @staticmethod
    def __generate_merge_polarity_isomorphic(motif_id: Union[int, str],
                                             polarities: list[list[str]],
                                             roles: list[tuple]) -> dict:
        isomorphic = []
        mappings = {}

        for curr_idx, pol in enumerate(polarities):
            sub_graph = _get_polarity_sub_graph(roles, pol)

            found = False
            found_idx = -1
            for found_idx, reps_graph in isomorphic:
                if isinstance(motif_id, str):
                    if _is_polarity_sim_isomorphic(reps_graph, sub_graph):
                        found = True
                        break
                else:
                    if nx.is_isomorphic(reps_graph,
                                        sub_graph,
                                        edge_match=lambda e1, e2: e1['polarity'] == e2['polarity']):
                        found = True
                        break

            curr_ = str(polarities[curr_idx])
            if found:
                mappings[curr_] = polarities[found_idx]

            if not found:
                mappings[curr_] = polarities[curr_idx]
                isomorphic.append((curr_idx, sub_graph))

        return mappings

    def merge_polarity_isomorphic_sub_graphs(
            self,
            motif_id: Union[int, str],
            polarity_frequencies: list[PolarityFrequencies]) -> list[PolarityFrequencies]:
        """
        merge polarity isomorphic sub graphs of the same base motif.
        e.g.: fan out "+ -" and "- +" are isomorphic, thus will be merged to a single pol motif
        and the other would be deleted from the returned list.
        :param motif_id: the original motif id (i.e.: not the polarity motif id)
        :param polarity_frequencies: a list of object containing frequencies and polarity list
        :return: merged polarity_frequencies list
        """
        del_idx_list = []
        pol_motif_mapping = self.polarity_mapping[motif_id]
        for curr_idx, pol_obj in enumerate(polarity_frequencies):
            iso_polarity = pol_motif_mapping[str(pol_obj.polarity)]
            if iso_polarity == pol_obj.polarity:
                continue

            merge_to_pol_obj = next(obj for obj in polarity_frequencies if obj.polarity == iso_polarity)
            merge_to_pol_obj.frequency += pol_obj.frequency
            merge_to_pol_obj.sub_graphs.extend(pol_obj.sub_graphs)
            del_idx_list.append(curr_idx)

        polarity_frequencies = [polarity_frequencies[i] for i in range(len(polarity_frequencies)) if
                                i not in del_idx_list]
        return polarity_frequencies

    def generate_polarity_isomorphic_k_sub_graphs(self):
        polarity_isomorphic_mapping = {}
        # TODO: extend with SIM keys
        for sub_id in tqdm(self.isomorphic_graphs.keys()):
            polarity_isomorphic_mapping[sub_id] = {}
            motif = create_base_motif(sub_id=sub_id, k=self.k)

            edges = len(motif.role_pattern)
            products = list(product(self.polarity_options, repeat=edges))
            polarities = [list(pol) for pol in products]
            polarity_isomorphic_mapping[sub_id] = self.__generate_merge_polarity_isomorphic(sub_id,
                                                                                            polarities,
                                                                                            motif.role_pattern)
        file_path = self.__get_polarity_isomorphic_k_file_name()
        self.__export_pol_iso_mapping(file_path, PolarityIsoMappingBinFile(mapping=polarity_isomorphic_mapping))

    @staticmethod
    def __export_iso_mapping(file_path: str, data: IsomorphicMappingBinFile):
        with open(file_path, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def __import_iso_mapping(file_path: str) -> IsomorphicMappingBinFile:
        with open(file_path, 'rb') as f:
            return IsomorphicMappingBinFile(**pickle.load(f))

    @staticmethod
    def __export_pol_iso_mapping(file_path: str, data: PolarityIsoMappingBinFile):
        with open(file_path, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def __import_pol_iso_mapping(file_path: str) -> PolarityIsoMappingBinFile:
        with open(file_path, 'rb') as f:
            return PolarityIsoMappingBinFile(**pickle.load(f))

    def __get_isomorphic_k_file_name(self):
        if self.allow_self_loops:
            file_path = f'{self.base_path}/k{self.k}_w_self_loops.bin'
        else:
            file_path = f'{self.base_path}/k{self.k}.bin'

        return file_path

    def __get_polarity_isomorphic_k_file_name(self):
        pol_name = 'cmpx_pol' if 'complex' in self.polarity_options else 'pol'
        if self.allow_self_loops:
            file_path = f'{self.base_path}/{pol_name}_k{self.k}_w_self_loops.bin'
        else:
            file_path = f'{self.base_path}/{pol_name}_k{self.k}.bin'
        return file_path


if __name__ == "__main__":
    iso_matcher = IsomorphicMotifMatch(k=3, polarity_options=['+', '-', 'complex'], allow_self_loops=False)
    # iso_matcher.generate_isomorphic_k_sub_graphs()
    # iso_matcher.generate_polarity_isomorphic_k_sub_graphs()
