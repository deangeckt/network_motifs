# https://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/uploads/SysBioCourse2018/exercise_2.pdf
from circuits import Circuits
from network import Network
from utils import Verbose


def analyze_network(file_path: str, name: str):
    print(f'Analyze {name} network:')
    with open(file_path, "r") as f:
        network = Network()
        network.load_graph(f.readlines())
        Circuits(network, verbose=Verbose.info).count_circuits()
    print()


if __name__ == "__main__":
    analyze_network("networks/toy_network.txt", "toy")
    analyze_network("networks/paper_exmp_network.txt", "paper example")
    analyze_network("networks/systems_biology_ex_network.txt", "uri alon course homework")
