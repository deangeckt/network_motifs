from enum import Enum
from typing import Optional, Union
from pydantic import BaseModel


class SubGraphAlgoName(str, Enum):
    specific = 'specific'
    mfinder_induced = 'mfinder_i'
    mfinder_none_induced = 'mfinder_ni'


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


class Motif(BaseModel):
    name: MotifName  # motif name
    id: int  # motif if / sub graph id
    adj_mat: str  # adjacency matrix

    # Optional: the motif is found in the real network:
    n_real: Optional[int] = 0  # number of appearances in the real network
    motif_criteria: Optional[MotifCriteriaResults] = None
    random_network_samples: Optional[list[int]] = []  # number of appearances of this motif id in the random networks
    sub_graphs: Optional[list[tuple]] = []  # all the isomorphic sub graphs appearances - in a tuple-edge format

    # a sorted dict of the nodes that appear in this motif
    # the key is either a neuron name or node id
    node_appearances: Optional[dict[Union[int, str], int]] = {}
