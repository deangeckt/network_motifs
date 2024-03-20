from circuits import Circuits
from motif_criteria import MotifCriteria
from network import Network
from network_randomizer import NetworkRandomizer
from simple_logger import Logger, LogLvl

logger = Logger(LogLvl.info)


def analyze_network(file_path: str, name: str):
    logger.info(f'Analyze {name} network:\n')
    with open(file_path, "r") as f:
        network = Network()
        network.load_from_txt(f.readlines())
        network.log_properties()
        # network.plot()
        network_circuits = Circuits(network).count_circuits()

        randomizer = NetworkRandomizer(network)
        random_network_amount = 1000

        logger.info(f"\nsearching network motifs using {random_network_amount} random networks")
        random_networks = randomizer.generate(amount=random_network_amount)
        random_network_circuits = [Circuits(network, use_logger=False).count_circuits() for network in random_networks]
        logger.toggle(True)

        for circuit in network_circuits:
            logger.info(f"\ncircuit {circuit}")
            random_network_samples = [rand_network[circuit] for rand_network in random_network_circuits]
            MotifCriteria().is_motif(network_circuits[circuit], random_network_samples)

    logger.info('\n')


if __name__ == "__main__":
    # analyze_network("networks/toy_network.txt", "toy")
    analyze_network("networks/paper_exmp_network.txt", "paper example")
    # analyze_network("networks/systems_biology_ex_network.txt", "uri alon course homework")
