import numpy as np
import scipy

from circuits import Circuits
from network import Network
from network_randomizer import NetworkRandomizer
from simple_logger import Logger, LogLvl

logger = Logger(LogLvl.debug)


def is_statistical_significant(n_real: int, random_network_samples: list[int], alpha=0.01) -> bool:
    n_rand = np.mean(random_network_samples)
    std = np.std(random_network_samples)
    if not std:
        logger.info('std is 0, cannot calculate z_score')
        return False

    z_score = (n_real - n_rand) / std
    p_value = scipy.stats.norm.sf(abs(z_score))
    logger.info(f'n_real: {n_real} n_rand: {n_rand} std: {std} z_score: {z_score}')
    is_significant = p_value < alpha
    logger.info(f'significant test with alpha: {alpha}: {is_significant}')
    return is_significant


def analyze_network(file_path: str, name: str):
    logger.info(f'Analyze {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_graph(f.readlines())
        network.log_properties()
        network_circuits = Circuits(network).count_circuits()

        randomizer = NetworkRandomizer(network)
        random_network_amount = 1000

        logger.info(f"\nsearching network motifs using {random_network_amount} random networks")
        random_networks = randomizer.generate(amount=random_network_amount)
        random_network_circuits = [Circuits(network, use_logger=False).count_circuits() for network in random_networks]
        logger.toggle(True)

        for circuit in network_circuits:
            logger.info(f"\n{circuit}:")
            random_network_samples = [rand_network[circuit] for rand_network in random_network_circuits]
            is_statistical_significant(network_circuits[circuit], random_network_samples)

    logger.info('\n')


if __name__ == "__main__":
    analyze_network("networks/toy_network.txt", "toy")
    # analyze_network("networks/paper_exmp_network.txt", "paper example")
    # analyze_network("networks/systems_biology_ex_network.txt", "uri alon course homework")
