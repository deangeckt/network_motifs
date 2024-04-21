from enum import Enum
from typing import Optional, Union

import numpy as np
from pydantic import BaseModel


class SubGraphAlgoName(str, Enum):
    specific = 'specific'
    mfinder_induced = 'mfinder_i'
    mfinder_none_induced = 'mfinder_ni'
    fanmod_esu = 'fanmod'


class RandomGeneratorAlgoName(str, Enum):
    markov_chain_switching = 'markov_chain'
    erdos_renyi = 'er'


class MotifType(str, Enum):
    motif = 'motif'
    anti_motif = 'anti-motif'
    none = 'none'


class MotifName(str, Enum):
    self_loop = 'self_loop'
    mutual_regulation = 'mutual_regulation'
    fan_out = 'fan_out'
    fan_in = 'fan_in'
    cascade = 'cascade'
    feed_forward = 'feed_forward'
    bi_fan = 'bi_fan'
    bi_parallel = 'bi_parallel'
    na = 'n/a'


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


class Motif(BaseModel):
    name: MotifName  # motif name
    id: int  # motif if / sub graph id
    adj_mat: np.ndarray  # adjacency matrix
    role_pattern: list[tuple]  # the adjacency matrix in an edge tuple format: e.g.: (a,b), (a,c)
    n_real: Optional[int] = 0  # number of appearances in the real network
    motif_criteria: Optional[MotifCriteriaResults] = None
    random_network_samples: Optional[list[int]] = []  # number of appearances of this motif id in the random networks
    sub_graphs: Optional[list[tuple[tuple]]] = []  # all the isomorphic sub graphs appearances - in a tuple-edge format
    node_roles: Optional[list[tuple]] = []  # list of (role, node) of all nodes participating in all the sub graphs

    # a sorted dict of the nodes that appear in this motif
    # the key is either a neuron name or node id
    node_appearances: Optional[dict[Union[int, str], int]] = {}

    polarity_frequencies: Optional[list[PolarityFrequencies]] = []

    class Config:
        arbitrary_types_allowed = True


class SubGraphSearchResult(BaseModel):
    # frequent sub graph list: key=motif id, value is the frequency
    fsl: dict[int, int]
    # same fsl, the value is the list of sub graphs
    fsl_fully_mapped: dict[int, list[tuple]]
