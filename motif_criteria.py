import numpy as np
import scipy

from utils.config import Config
from utils.simple_logger import Logger


class MotifCriteria:
    def __init__(self):
        self.logger = Logger()
        config = Config()

        self.alpha = float(config.get_property('motif_criteria', 'alpha'))
        self.n_real_min = int(config.get_property('motif_criteria', 'n_real_min'))
        self.larger_than_rand_factor = float(config.get_property('motif_criteria', 'larger_than_rand_factor'))

    def __is_statistical_significant_rand_network(self, n_real: int, random_network_samples: list[int]) -> bool:
        n_rand = np.mean(random_network_samples)
        std = np.std(random_network_samples)
        if not std:
            self.logger.info('std is 0, cannot calculate z_score')
            return False

        z_score = (n_real - n_rand) / std
        p_value = scipy.stats.norm.sf(abs(z_score))
        self.logger.info(f'n_real: {n_real} n_rand: {n_rand} std: {round(std, 2)} z_score: {round(z_score, 2)}')
        is_significant = p_value < self.alpha
        self.logger.info(f'significant test; p<alpha: {round(p_value, 2)} < {self.alpha}: {is_significant}')
        return is_significant

    def __is_larger_than_rand(self, n_real: int, random_network_samples: list[int]) -> bool:
        n_rand = np.mean(random_network_samples)
        is_larger = (n_real - n_rand) > (self.larger_than_rand_factor * n_rand)
        self.logger.info(
            f'is larger than random test; real-rand>0.1*rand: '
            f'{round(n_real - n_rand, 2)} > {round(self.larger_than_rand_factor * n_rand, 2)}: {is_larger}')

        return is_larger

    def is_motif(self, n_real: int, random_network_samples: list[int]) -> bool:
        self.logger.info('checking motif criteria:')
        is_significant = self.__is_statistical_significant_rand_network(n_real, random_network_samples)
        is_larger_than_rand = self.__is_larger_than_rand(n_real, random_network_samples)

        is_min = n_real >= self.n_real_min
        self.logger.info(f'is min test: {n_real} >= {self.n_real_min}: {is_min}')

        is_motif = is_min & is_significant & is_larger_than_rand
        self.logger.info(f'motif criteria: {is_motif}')
        return is_motif
