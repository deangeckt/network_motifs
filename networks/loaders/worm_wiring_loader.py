import numpy as np
import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from networks.network import Network


class WormWiringLoader(NetworkLoaderStrategy):
    def __init__(self, args):
        super().__init__(args)
        self.amount = 272

    def load(self, *args) -> Network:
        """
        https://wormwiring.org/pages/adjacency.html xlsx format
        Cook, (2019). Whole-animal connectomes of both Caenorhabditis elegans sexes.
        """
        xlsx_path, sheet_name = args
        amount = self.amount

        row_start = 60 if 'SI 2' in xlsx_path else 23  # SI 5
        col_start = 60 if 'SI 2' in xlsx_path else 53  # SI 5

        xls = pd.ExcelFile(xlsx_path)
        df = xls.parse(sheet_name, header=None)

        names_idx = 2
        adj_mat = np.array(df.iloc[row_start:row_start + amount, col_start:col_start + amount])
        neurons_col = list(df.iloc[row_start:row_start + amount, names_idx])
        neurons_row = list(df.iloc[names_idx, col_start:col_start + amount])

        if not neurons_row == neurons_col:
            raise Exception(f'invalid xlsx file or indices')

        self.neuron_names = neurons_col

        for i in range(len(adj_mat)):
            for j in range(len(adj_mat)):
                synapses = adj_mat[i, j]
                if np.isnan(synapses) or synapses == 0:
                    continue
                self._load_synapse(i, j, num_of_synapse=synapses, polarity=None)

        return self._copy_network_params()
