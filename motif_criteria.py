import numpy as np
import scipy

from utils.config import Config
from utils.simple_logger import Logger

from enum import Enum


class MotifType(str, Enum):
    motif = 'motif'
    anti_motif = 'anti-motif'
    none = 'none'


class MotifCriteria:
    def __init__(self):
        self.logger = Logger()
        config = Config()

        self.alpha = float(config.get_property('motif_criteria', 'alpha'))

        # The uniqueness value is the number of times a subgraph appears in the
        # network with completely disjoint groups of nodes
        self.uniqueness_threshold = int(config.get_property('motif_criteria', 'uniqueness_threshold'))

        # Mfactor in the original paper
        self.frequency_threshold = float(config.get_property('motif_criteria', 'frequency_threshold'))

    def __is_statistically_significant(self, n_real: int, random_network_samples: list[int]) -> bool:
        n_rand = np.mean(random_network_samples)
        self.logger.info(f'n_real: {n_real} n_rand: {round(n_rand, 2)}')

        std = np.std(random_network_samples)
        if not std:
            self.logger.info('std is 0, cannot calculate z_score')
            return False

        z_score = (n_real - n_rand) / std
        p_value = scipy.stats.norm.sf(abs(z_score))
        self.logger.info(f'std: {round(std, 2)} z_score: {round(z_score, 2)}')

        is_significant = p_value < self.alpha
        self.logger.info(f'\tsignificant test; p<alpha: {round(p_value, 2)} < {self.alpha}: {is_significant}')

        return is_significant

    def __is_frequency(self, n_real: int, random_network_samples: list[int]) -> bool:
        n_rand = np.mean(random_network_samples)
        is_freq = (n_real - n_rand) > (self.frequency_threshold * n_rand)
        self.logger.info(
            f'\tfrequency threshold test; real-rand > {self.frequency_threshold}*rand: '
            f'{round(n_real - n_rand, 2)} > {round(self.frequency_threshold * n_rand, 2)}: {is_freq}')

        return is_freq

    def __is_anti_frequency(self, n_real: int, random_network_samples: list[int]) -> bool:
        n_rand = np.mean(random_network_samples)
        is_freq = (n_rand - n_real) > (self.frequency_threshold * n_rand)
        self.logger.info(
            f'\tanti frequency threshold test; rand-real > {self.frequency_threshold}*rand: '
            f'{round(n_rand - n_real, 2)} > {round(self.frequency_threshold * n_rand, 2)}: {is_freq}')

        return is_freq

    def is_motif(self, n_real: int, random_network_samples: list[int]) -> MotifType:
        is_significant = self.__is_statistically_significant(n_real, random_network_samples)
        is_freq = self.__is_frequency(n_real, random_network_samples)

        # TODO: fix uniqueness test
        is_larger_than_min = n_real >= self.uniqueness_threshold
        self.logger.info(f'\tis uniqueness test; {n_real} >= {self.uniqueness_threshold}: {is_larger_than_min}')

        is_motif = is_larger_than_min & is_significant & is_freq
        res = MotifType.motif if is_motif else MotifType.none

        # check for anti motif
        if not is_motif:
            if self.__is_anti_frequency(n_real, random_network_samples) and is_significant:
                res = MotifType.anti_motif

        self.logger.info(f'Motif criteria result: {res}')
        return res
