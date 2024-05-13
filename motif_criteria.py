import numpy as np
import scipy

from subgraphs.sub_graphs_utils import get_number_of_disjoint_group_nodes
from utils.logs import log_motif_criteria_args
from utils.simple_logger import Logger
from utils.types import MotifCriteriaResults, MotifType, Motif, MotifCriteriaArgs


class MotifCriteria:
    def __init__(self, args: MotifCriteriaArgs):
        self.logger = Logger()

        self.alpha = args.alpha

        # The uniqueness value is the number of times a subgraph appears in the
        # network with completely disjoint groups of nodes
        self.uniqueness_threshold = args.uniqueness_threshold
        self.use_uniqueness = args.use_uniq_criteria

        # Mfactor in the original paper
        self.frequency_threshold = args.frequency_threshold

        log_motif_criteria_args(args)

    def __calculate_statistics(self, n_real: int,
                               random_network_samples: list[int]) -> MotifCriteriaResults:

        n_rand = np.mean(random_network_samples)
        std = np.std(random_network_samples)
        if not std:
            self.logger.debug('\tstd is 0, cannot calculate z_score')
            return MotifCriteriaResults(
                n_rand=n_rand,
                n_real=n_real,
                is_statistically_significant=False
            )

        z_score = (n_real - n_rand) / std
        p_value = scipy.stats.norm.sf(abs(z_score))
        is_significant = p_value < self.alpha

        self.logger.debug(f'\tstd: {round(std, 2)} z_score: {round(z_score, 2)}')
        self.logger.debug(f'\tsignificant test; p<alpha: {round(p_value, 2)} < {self.alpha}: {is_significant}')

        return MotifCriteriaResults(
            z_score=z_score,
            std=std,
            p_value=p_value,
            n_real=n_real,
            n_rand=n_rand,
            is_statistically_significant=is_significant
        )

    def __is_frequency(self, partial_res: MotifCriteriaResults):
        n_rand = partial_res.n_rand
        n_real = partial_res.n_real

        is_freq = (n_real - n_rand) > (self.frequency_threshold * n_rand)
        self.logger.debug(
            f'\tfrequency threshold test; real-rand > {self.frequency_threshold}*rand: '
            f'{round(n_real - n_rand, 2)} > {round(self.frequency_threshold * n_rand, 2)}: {is_freq}')

        partial_res.is_motif_frequent = is_freq

    def __is_anti_frequency(self, partial_res: MotifCriteriaResults):
        n_rand = partial_res.n_rand
        n_real = partial_res.n_real

        is_freq = (n_rand - n_real) > (self.frequency_threshold * n_rand)
        self.logger.debug(
            f'\tanti frequency threshold test; rand-real > {self.frequency_threshold}*rand: '
            f'{round(n_rand - n_real, 2)} > {round(self.frequency_threshold * n_rand, 2)}: {is_freq}')

        partial_res.is_anti_motif_frequent = is_freq

    def __is_uniq(self, partial_res: MotifCriteriaResults, sub_graphs: list):
        """
        :param sub_graphs: a list of the motif candidate sub graphs edges
        :return: if configured to use the uniqueness test
        check if the motif candidate uniq amount of sub graphs >= uniqueness_threshold
        """
        uniq_n_real = get_number_of_disjoint_group_nodes(sub_graphs)
        is_uniq = uniq_n_real >= self.uniqueness_threshold
        self.logger.debug(f'\tis uniqueness test; {uniq_n_real} >= {self.uniqueness_threshold}: {is_uniq}')

        partial_res.is_uniq = is_uniq if self.use_uniqueness else 'n/a'
        partial_res.uniq = uniq_n_real

    def is_motif(self, motif: Motif) -> MotifCriteriaResults:
        self.logger.debug(f'Motif criteria for {motif.id}')
        motif_results = self.__calculate_statistics(motif.n_real, motif.random_network_samples)

        self.__is_frequency(motif_results)
        self.__is_anti_frequency(motif_results)
        self.__is_uniq(motif_results, motif.sub_graphs)

        is_motif: bool = (motif_results.is_uniq and
                          motif_results.is_statistically_significant and
                          motif_results.is_motif_frequent)

        motif_results.is_motif = MotifType.motif if is_motif else MotifType.none

        # check for anti motif
        if not is_motif:
            if motif_results.is_anti_motif_frequent and motif_results.is_statistically_significant:
                motif_results.is_motif = MotifType.anti_motif

        return motif_results
