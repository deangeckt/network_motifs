from argparse import Namespace
from enum import Enum
from typing import Optional, Union, TypedDict

import numpy as np
from pydantic import BaseModel


class NetworkInputType(str, Enum):
    simple_adj_txt = 'simple_adj_txt'
    worm_wiring_xlsx = 'worm_wiring_xlsx'
    polarity_xlsx = 'polarity_xlsx'
    durbin_txt = 'durbin_txt'
    graph = 'graph'


class NetworkLoaderArgs(BaseModel):
    synapse_threshold: int
    filter_polarity: Optional[list[str]] = ['+', '-']
    filter_prim_nt: Optional[list[Union[str, int]]] = ['GABA', 'Glu', 'ACh', 0]
    filter_syn_type: Optional[str] = 'chem'
    filter_sex_type: Optional[str] = 'herm'
    filter_nerve_ring_neurons: Optional[bool] = False


class SubGraphAlgoName(str, Enum):
    specific = 'specific'
    mfinder_induced = 'mfinder_i'
    mfinder_none_induced = 'mfinder_ni'
    fanmod_esu = 'fanmod'
    triadic_census = 'triadic_census'


class RandomGeneratorAlgoName(str, Enum):
    markov_chain_switching = 'markov_chain'
    nerve_ring_markov_chain_switching = 'nerve_ring_markov_chain'
    erdos_renyi = 'erdos_renyi'
    barabasi = 'barabasi'


class MotifType(str, Enum):
    motif = 'motif'
    anti_motif = 'anti-motif'
    none = 'none'


class MotifName(str, Enum):
    self_loop = 'self_loop'
    mutual_regulation = 'mutual_regulation'
    fan_out = 'fan_out'  # a.k.a sim_2
    fan_in = 'fan_in'
    cascade = 'cascade'
    feed_forward = 'feed_forward'
    bi_fan = 'bi_fan'
    bi_parallel = 'bi_parallel'
    sim_3 = 'sim_3'
    na = 'n/a'


class MotifCriteriaArgs(BaseModel):
    alpha: float
    uniqueness_threshold: int
    use_uniq_criteria: bool
    frequency_threshold: float


class MotifCriteriaResults(BaseModel):
    n_real: int
    is_statistically_significant: Optional[bool] = False
    n_rand: Optional[float] = 0.0
    z_score: Optional[float] = 0.0
    std: Optional[float] = 0.0
    p_value: Optional[float] = 0.0
    uniq: Optional[int] = -1
    is_motif_frequent: Optional[bool] = False
    is_anti_motif_frequent: Optional[bool] = False
    is_uniq: Optional[Union[bool, str]] = False
    is_motif: Optional[MotifType] = MotifType.none


class PolarityFrequencies(BaseModel):
    frequency: int
    polarity: list[str]
    sub_graphs: list


class Motif(BaseModel):
    name: MotifName  # motif name
    id: Union[int, str]  # motif if / sub graph id
    adj_mat: np.ndarray  # adjacency matrix
    role_pattern: list[tuple]  # the adjacency matrix in an edge tuple format: e.g.: (a,b), (a,c)
    n_real: Optional[int] = 0  # number of appearances in the real network
    motif_criteria: Optional[MotifCriteriaResults] = None
    random_network_samples: Optional[list[int]] = []  # number of appearances of this motif id in the random networks

    # all the isomorphic sub graphs appearances - in a tuple-edge format: (s,t,polarity)
    sub_graphs: Optional[list[tuple[tuple]]] = []
    # dict of dicts, key=role. value = dict where keys are node name and value are their freq.
    node_roles: Optional[dict] = {}

    # a sorted dict of the nodes that appear in this motif
    # the key is either a neuron name or node id
    node_appearances: Optional[dict[Union[int, str], int]] = {}

    polarity_motifs: Optional[list['Motif']] = []
    polarity: Optional[list[str]] = []

    class Config:
        arbitrary_types_allowed = True


class SubGraphSearchResult(BaseModel):
    # frequent sub graph list: key=motif id, value is the frequency
    fsl: dict[Union[str, int], int]
    # same fsl, the value is the list of sub graphs
    fsl_fully_mapped: dict[Union[str, int], list[tuple]]


class LargeSubGraphSearchResult(SubGraphSearchResult):
    adj_mat: dict[str, np.ndarray]

    class Config:
        arbitrary_types_allowed = True


class SearchResultBinaryFile(TypedDict):
    args: Namespace
    motifs: dict[Union[str, int], Motif]
