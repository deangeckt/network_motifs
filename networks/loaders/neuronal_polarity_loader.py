import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from networks.network import Network
from utils.config import Config


class NeuronalPolarityLoader(NetworkLoaderStrategy):
    def __init__(self):
        """
        xlsx files from the paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P.
        Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical
        synapse neuronal connectome network
        """
        super().__init__()
        config = Config()

        # polarity configuration
        self.src_col = 0
        self.tar_col = 3
        self.weight_col = 4
        self.edge_type_col = 5
        self.polarity_col = 16
        self.prim_nt_col = 1

        # polarity options [+, -, no pred, complex]
        self.filter_polarity = config.get_string_list('polarity', 'filter_polarity')
        # primary neurotransmitter options [GABA, Glu, ACh, 0] (0 is an int)
        self.filter_prim_nt = config.get_string_list('polarity', 'filter_prim_nt')

    def load(self, *args) -> Network:
        xlsx_path, sheet_name = args
        self.use_polarity = True
        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(sheet_name, header=None)

        # filter
        self.logger.info(f'\nFiltering Neurons with polarity: {self.filter_polarity}')
        df = df[df[self.polarity_col].isin(self.filter_polarity)]
        self.logger.info(f'Filtering Neurons with primary neurotransmitter: {self.filter_prim_nt}\n')
        df = df[df[self.prim_nt_col].isin(self.filter_prim_nt)]

        src_neurons_names = df.iloc[:, self.src_col]
        tar_neurons_names = df.iloc[:, self.tar_col]
        edge_weights = df.iloc[:, self.weight_col]  # these are the amount of synapses
        polarity = df.iloc[:, self.polarity_col]

        self.neuron_names = list(set(src_neurons_names) | set(tar_neurons_names))
        neurons_indices = {ss: i for i, ss in enumerate(self.neuron_names)}

        for v1, v2, w, p in zip(src_neurons_names, tar_neurons_names, edge_weights, polarity):
            polarity_edge = 1 if p == '+' else -1
            self._load_synapse(neurons_indices[v1], neurons_indices[v2], w, polarity_edge)

        return self._copy_network_params()
