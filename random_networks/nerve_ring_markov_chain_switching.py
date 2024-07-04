import random

import pandas as pd
from networkx import DiGraph

from networks.network import Network
from random_networks.markov_chain_switching import MarkovChainSwitching
from collections import defaultdict

nerve_ring_neurons = [
    'ADAR',
    'ADAL',
    'ADEL',
    'ADER',
    'ADFL',
    'ADFR',
    'ADLL',
    'ADLR',
    'AFDL',
    'AFDR',
    'AIAL',
    'AIAR',
    'AIBL',
    'AIBR',
    'AIML',
    'AIMR',
    'AINL',
    'AINR',
    'AIYL',
    'AIYR',
    'AIZL',
    'AIZR',
    'ALA',
    'ALML',
    'ALMR',
    'ALNL',
    'ALNR',
    'AQR',
    'ASEL',
    'ASER',
    'ASGL',
    'ASGR',
    'ASHL',
    'ASHR',
    'ASIL',
    'ASIR',
    'ASJL',
    'ASJR',
    'ASKL',
    'ASKR',
    'AUAL',
    'AUAR',
    'AVAL',
    'AVAR',
    'AVBL',
    'AVBR',
    'AVDL',
    'AVDR',
    'AVEL',
    'AVER',
    'AVFL',
    'AVFR',
    'AVHL',
    'AVHR',
    'AVJL',
    'AVJR',
    'AVKL',
    'AVKR',
    'AVL',
    'AVM',
    'AWAL',
    'AWAR',
    'AWBL',
    'AWBR',
    'AWCL',
    'AWCR',
    'BAGL',
    'BAGR',
    'BDUL',
    'BDUR',
    'CEPDL',
    'CEPDR',
    'CEPVL',
    'CEPVR',
    'DVA',
    'DVC',
    'FLPL',
    'FLPR',
    'HSNL',
    'HSNR',
    'IL1DL',
    'IL1DR',
    'IL1L',
    'IL1R',
    'IL1VL',
    'IL1VR',
    'IL2DL',
    'IL2DR',
    'IL2L',
    'IL2R',
    'IL2VL',
    'IL2VR',
    'OLLL',
    'OLLR',
    'OLQDL',
    'OLQDR',
    'OLQVL',
    'OLQVR',
    'PLNL',
    'PLNR',
    'PVCL',
    'PVCR',
    'PVNL',
    'PVNR',
    'PVPL',
    'PVPR',
    'PVQL',
    'PVQR',
    'PVR',
    'PVT',
    'RIAL',
    'RIAR',
    'RIBL',
    'RIBR',
    'RICL',
    'RICR',
    'RID',
    'RIFL',
    'RIFR',
    'RIGL',
    'RIGR',
    'RIH',
    'RIML',
    'RIMR',
    'RIPL',
    'RIPR',
    'RIR',
    'RIS',
    'RIVL',
    'RIVR',
    'RMDDL',
    'RMDDR',
    'RMDL',
    'RMDR',
    'RMDVL',
    'RMDVR',
    'RMED',
    'RMEL',
    'RMER',
    'RMEV',
    'RMFL',
    'RMFR',
    'RMGL',
    'RMGR',
    'RMHL',
    'RMHR',
    'SAADL',
    'SAADR',
    'SAAVL',
    'SAAVR',
    'SDQL',
    'SDQR',
    'SIADL',
    'SIADR',
    'SIAVL',
    'SIAVR',
    'SIBDL',
    'SIBDR',
    'SIBVL',
    'SIBVR',
    'SMBDL',
    'SMBDR',
    'SMBVL',
    'SMBVR',
    'SMDDL',
    'SMDDR',
    'SMDVL',
    'SMDVR',
    'URADL',
    'URADR',
    'URAVL',
    'URAVR',
    'URBL',
    'URBR',
    'URXL',
    'URXR',
    'URYDL',
    'URYDR',
    'URYVL',
    'URYVR',
    'VB01'
]


class NerveRingMarkovChainSwitching(MarkovChainSwitching):
    """
    Nerve ring distances are based on the matlab file from:
    Zaslaver et al., 2022:
    "The synaptic organization in the Caenorhabditis elegans
    neural network suggests significant local compartmentalized
    computations"
    """

    def init_nerve_ring_allow_list(self):
        xls = pd.ExcelFile('random_networks/nerve_ring_distances.xlsx')
        df = xls.parse('distance', header=None)

        for i, neuron in enumerate(nerve_ring_neurons):
            for j in range(len(nerve_ring_neurons)):
                if i == j:
                    continue
                dist = df.iloc[i][j]
                if dist < self.DISTANCE_TH:
                    self.nerve_ring_allow[neuron].add(nerve_ring_neurons[j])

    def __init__(self, network: Network, switch_factor: int):
        super().__init__(network, switch_factor)

        self.DISTANCE_TH = 0.5
        self.neuron_names = self.network.neuron_names
        self.nerve_ring_allow = defaultdict(set)
        self.init_nerve_ring_allow_list()

    def nerve_ring_constrain(self, x1: str, y1: str, x2: str, y2: str) -> bool:
        # in case x is not in the nerve ring (not in the dict), we allow the switch.
        if y2 not in self.nerve_ring_allow.get(x1, {y2}):
            return False

        if y1 not in self.nerve_ring_allow.get(x2, {y1}):
            return False

        return True

    def _markov_chain(self) -> DiGraph:
        graph: DiGraph = self.network.graph.copy()

        for _ in range(self.markov_chain_num_iterations):
            e1, e2 = random.sample(graph.edges, 2)
            x1, y1 = e1
            x2, y2 = e2

            if not self.switching_constrain(graph, x1, x2, y1, y2):
                continue

            if not self.nerve_ring_constrain(self.neuron_names[x1], self.neuron_names[x2],
                                             self.neuron_names[y1], self.neuron_names[y2]):
                continue

            self.success_switch += 1
            self.switch_foo(graph, x1, y1, x2, y2)

        return graph
