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
    sum_ = 0
    warn_amount_max = iterations / 2
    warn_sum_ = 0
    while sum_ < iterations:
        samples = np.random.choice(motif.random_network_samples, size=iteration_sample, replace=True)
        n_rand = np.mean(samples)
        std = np.std(samples)
        if not std:
            warn_sum_ += 1
            if warn_sum_ == warn_amount_max:
                print(f'motif id: {motif.id} err: std was 0, {warn_amount_max} times! z scores len is: {len(z_scores)}')
                return z_scores
            continue
        z_scores.append((motif.n_real - n_rand) / std)
        sum_ += 1
    return z_scores
