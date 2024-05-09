import numpy as np

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from networks.network import Network


class DurbinFileLoader(NetworkLoaderStrategy):
    def __init__(self, args):
        """
        https://www.wormatlas.org/neuronalwiring.html - Neuronal Connectivity I: by R. Durbin 1986
        : param filter_syn_type: either 'chem', 'gap', 'all
        : param filter_recon: either: one of 'JSH', 'N2U'
        """
        super().__init__(args)

        # 'chem', 'gap', 'all
        self.filter_syn_type = args.durbin_filter_syn_type

        # 'JSH: L4 male', 'N2U: hermaphrodite adult'
        self.filter_recon = args.durbin_filter_recon

        self.logger.info(f'Filtering Neurons of: {self.filter_recon}')
        self.logger.info(f'Filtering Synapses of type: {self.filter_syn_type}')

    @staticmethod
    def __append_adj_mat(adj_mat: np.ndarray, v1: int, v2: int, num_of_synapses: int):
        if not np.isnan(adj_mat[v1, v2]):
            adj_mat[v1, v2] += num_of_synapses
        else:
            adj_mat[v1, v2] = num_of_synapses

    def load(self, *args) -> Network:
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

        self.neuron_names = list(neurons_names)
        neurons = {n: i for i, n in enumerate(self.neuron_names)}
        N = len(neurons)
        adj_mat = np.empty((N, N))
        adj_mat.fill(np.nan)

        for n1, n2, synapse_type, _, num_of_synapses in data:
            if 'Receive' in synapse_type:
                self.__append_adj_mat(adj_mat=adj_mat, v1=neurons[n2], v2=neurons[n1],
                                      num_of_synapses=int(num_of_synapses))
            elif synapse_type == 'Gap_junction':
                self.__append_adj_mat(adj_mat=adj_mat, v1=neurons[n1], v2=neurons[n2],
                                      num_of_synapses=int(num_of_synapses))
                self.__append_adj_mat(adj_mat=adj_mat, v1=neurons[n2], v2=neurons[n1],
                                      num_of_synapses=int(num_of_synapses))
            else:
                self.__append_adj_mat(adj_mat=adj_mat, v1=neurons[n1], v2=neurons[n2],
                                      num_of_synapses=int(num_of_synapses))

        for i in range(len(adj_mat)):
            for j in range(len(adj_mat)):
                synapses = adj_mat[i, j]
                if np.isnan(synapses) or synapses == 0:
                    continue
                self._load_synapse(i, j, num_of_synapse=synapses, polarity=None)

        return self._copy_network_params()
