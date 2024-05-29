import numpy as np

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.types import NetworkLoaderArgs


class DurbinFileLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        https://www.wormatlas.org/neuronalwiring.html - Neuronal Connectivity I: by R. Durbin 1986
        """
        super().__init__(args)

        # 'chem', 'gap', 'all'
        self.filter_syn_type = args.filter_syn_type

        # 'JSH': L4 male', 'N2U: hermaphrodite adult'
        self.filter_recon = 'N2U' if args.filter_sex_type == 'herm' else 'JSH'

        self.logger.info(f'Filtering Neurons of: {self.filter_recon}')
        self.logger.info(f'Filtering Synapses of type: {self.filter_syn_type}')

    def load(self, *args):
        neurons_names = set()
        data = []
        file_path = args[0]

        with (open(file_path, "r") as f):
            for line in f.readlines():
                line = tuple(line.strip().split())
                if len(line) == 5:
                    n1, n2, synapse_type, reconstruction, num_of_synapses = line
                elif len(line) == 4:
                    n1, n2, synapse_type, num_of_synapses = line
                    reconstruction = 'none'
                else:
                    raise Exception('invalid durbin file')
                data.append((n1, n2, synapse_type, reconstruction, int(num_of_synapses)))

        data = list(filter(lambda x: x[3] == self.filter_recon, data))
        gap_data = list(filter(lambda x: x[2] == 'Gap_junction', data))
        chem_data = list(filter(lambda x: x[2] != 'Gap_junction', data))

        if self.filter_syn_type == 'all':
            data = gap_data + chem_data
        elif self.filter_syn_type == 'gap':
            data = gap_data
        elif self.filter_syn_type == 'chem':
            data = chem_data
        else:
            raise Exception(f'invalid syn type: {self.filter_syn_type}')

        for n1, n2, _, _, _ in data:
            neurons_names.add(n1)
            neurons_names.add(n2)

        self.network.neuron_names = list(neurons_names)
        neurons = {n: i for i, n in enumerate(self.network.neuron_names)}

        for n1, n2, synapse_type, _, num_of_synapses in data:
            if 'Receive' in synapse_type:
                self._append_edge(v1=neurons[n2], v2=neurons[n1], synapse=int(num_of_synapses), gap=0)
            elif synapse_type == 'Gap_junction':
                self._append_edge(v1=neurons[n1], v2=neurons[n2], synapse=0, gap=int(num_of_synapses))
                self._append_edge(v1=neurons[n2], v2=neurons[n1], synapse=0, gap=int(num_of_synapses))
            else:
                self._append_edge(v1=neurons[n1], v2=neurons[n2], synapse=int(num_of_synapses), gap=0)
