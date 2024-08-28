from utils.types import Motif
import numpy as np


def z_score_bootstrap(motif: Motif, iterations=10000, iteration_sample=1000) -> list[float]:
    """
    calculate z score distribution, given n_real and sampled random networks of a given motif
    :param motif: the motif to which we calculate the z score distribution
    :param iterations: number of bootstrap iterations; will result in a distribution of size 'iterations'
    :param iteration_sample: sampling n random networks with repeats
    :return: z score distribution (list of floats)
    """
    z_scores = []
    for _ in range(iterations):
        samples = np.random.choice(motif.random_network_samples, size=iteration_sample, replace=True)
        n_rand = np.mean(samples)
        std = np.std(samples)
        if not std:
            continue
        z_scores.append((motif.n_real - n_rand) / std)
    return z_scores
