import collections
import os
import pickle
from collections import defaultdict
from itertools import product
from typing import Union, TypedDict

import networkx as nx
from networkx.algorithms import isomorphism
from tqdm import tqdm

from large_subgraphs.single_input_moudle import get_sim_adj_mat
from utils.sub_graphs import get_sub_graph_from_id, create_base_motif, create_sim_motif
from utils.types import PolarityFrequencies, Motif


class IsomorphicMappingBinFile(TypedDict):
    isomorphic_mapping: dict
    isomorphic_graphs: dict


class PolarityIsoMappingBinFile(TypedDict):
    mapping: dict


def match_two_fsl_id_lists(src_fsl_ids: list[int], tar_fsl_ids: list[int], k: int) -> dict[int, int]:
    """
    Used when use_isomorphic_mapping flag is off, i.e.: larger K's.
    #TODO: does this support the polarity network with large k?
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
        self.k = k
        self.allow_self_loops = allow_self_loops
        self.polarity_options = polarity_options

        self.min_sim_ctrl_size = 3
        self.max_sim_ctrl_size = 12

        self.base_path = 'isomorphic/mapping'

        # Load iso mapping
        file_path = self.__get_isomorphic_k_file_name()
        if not os.path.isfile(file_path):
            print(f'isomorphic file for k={k} does not exist')
            return

        bin_file: IsomorphicMappingBinFile = self.__import_iso_mapping(file_path)
        self.isomorphic_mapping = bin_file['isomorphic_mapping']
        self.isomorphic_graphs = bin_file['isomorphic_graphs']

        # Load polarity iso mapping
        if not polarity_options:
            return

        file_path = self.__get_polarity_isomorphic_k_file_name()
        if not os.path.isfile(file_path):
            print(f'polarity isomorphic {self.polarity_options} file for k={k} does not exist')
            return
        bin_file: PolarityIsoMappingBinFile = self.__import_pol_iso_mapping(file_path)
        self.polarity_mapping = bin_file['mapping']

        # Load SIM polarity iso mapping
        file_path = self.__get_polarity_isomorphic_sim_file_name()
        if not os.path.isfile(file_path):
            print(f'sim polarity isomorphic {self.polarity_options} '
                  f'file for control={self.max_sim_ctrl_size} does not exist')
            return
        bin_file: PolarityIsoMappingBinFile = self.__import_pol_iso_mapping(file_path)
        self.sim_mapping = bin_file['mapping']

    def merge_polarity_isomorphic_frequencies(
            self,
            motif_id: Union[int, str],
            polarity_frequencies: list[PolarityFrequencies]) -> list[PolarityFrequencies]:
        """
        merge polarity isomorphic sub graphs frequencies of the same base motif.
        e.g.: fan out "+ -" and "- +" are isomorphic, thus will be merged to a single pol motif
        and the other would be deleted from the returned list.
        :param motif_id: the original motif id (i.e.: not the polarity motif id)
        :param polarity_frequencies: a list of object containing frequencies and polarity list
        :return: merged polarity_frequencies list
        """
        del_idx_list = []
        if isinstance(motif_id, str):
            pol_motif_mapping = self.sim_mapping[motif_id]
        else:
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

    def generate_isomorphic_k_sub_graphs(self):
        """
        export: isomorphic_mapping: a dict where each key is a motif / sub graph id and
        the value is the smallest motif id of the same isomorphic set of sub graphs.
        export isomorphic_graphs: a dict where each key is the smallest motif id and the value
        is the list of all the motif ids that are isomorphic to it.

        not scalable for k >= 5.
        """
        isomorphic = defaultdict(dict)

        possible_options = (2 ** (self.k ** 2))
        for sub_id in tqdm(range(possible_options)):
            sub_graph = get_sub_graph_from_id(decimal=sub_id, k=self.k)

            if not self.allow_self_loops and list(nx.selfloop_edges(sub_graph)):
                continue

            un_dir_sub_graph = nx.Graph(sub_graph)
            # remove the not connected cases, e.g.: no edges at all-sub graph
            if not nx.is_connected(un_dir_sub_graph):
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
    def __generate_and_merge_polarity_isomorphic(motif_id: Union[int, str],
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

    def __generate_polarity_isomorphic_for_motif(self, motif: Motif) -> dict:
        edges = len(motif.role_pattern)
        products = list(product(self.polarity_options, repeat=edges))
        polarities = [list(pol) for pol in products]
        return self.__generate_and_merge_polarity_isomorphic(motif.id, polarities, motif.role_pattern)

    def generate_polarity_isomorphic_k_sub_graphs(self):
        polarity_isomorphic_mapping = {}

        for sub_id in tqdm(self.isomorphic_graphs.keys()):
            polarity_isomorphic_mapping[sub_id] = {}
            motif = create_base_motif(sub_id=sub_id, k=self.k)
            polarity_isomorphic_mapping[sub_id] = self.__generate_polarity_isomorphic_for_motif(motif)

        file_path = self.__get_polarity_isomorphic_k_file_name()
        self.__export_pol_iso_mapping(file_path, PolarityIsoMappingBinFile(mapping=polarity_isomorphic_mapping))

    def generate_polarity_isomorphic_sim_graphs(self):
        polarity_isomorphic_mapping = {}

        for control_size in tqdm(range(self.min_sim_ctrl_size, self.max_sim_ctrl_size + 1)):
            sim_id = f'SIM_{control_size}'
            polarity_isomorphic_mapping[sim_id] = {}
            motif = create_sim_motif(sim_id=sim_id, adj_mat=get_sim_adj_mat(control_size))
            polarity_isomorphic_mapping[sim_id] = self.__generate_polarity_isomorphic_for_motif(motif)

        file_path = self.__get_polarity_isomorphic_sim_file_name()
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
            file_path = f'{self.base_path}/k{self.k}_{pol_name}_w_self_loops.bin'
        else:
            file_path = f'{self.base_path}/k{self.k}_{pol_name}.bin'
        return file_path

    def __get_polarity_isomorphic_sim_file_name(self):
        pol_name = 'cmpx_pol' if 'complex' in self.polarity_options else 'pol'
        file_path = f'{self.base_path}/sim{self.max_sim_ctrl_size}_{pol_name}.bin'
        return file_path


if __name__ == "__main__":
    # iso_matcher = IsomorphicMotifMatch(k=4, polarity_options=['+', '-'], allow_self_loops=False)
    iso_matcher = IsomorphicMotifMatch(k=3, polarity_options=[], allow_self_loops=True)

    iso_matcher.generate_isomorphic_k_sub_graphs()
    # iso_matcher.generate_polarity_isomorphic_k_sub_graphs()
    # iso_matcher.generate_polarity_isomorphic_sim_graphs()
