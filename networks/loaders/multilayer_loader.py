import pandas as pd

from networks.loaders.network_loader_strategy import NetworkLoaderStrategy
from utils.polarity_counter import count_network_polarity_ratio
from utils.types import NetworkLoaderArgs


class MultilayerConnectomeLoader(NetworkLoaderStrategy):
    def __init__(self, args: NetworkLoaderArgs):
        """
        paper: Bentley B, Branicky R, Barnes CL, Chew YL, Yemini E, et al. (2016)
        The Multilayer Connectome of Caenorhabditis elegans
        """
        super().__init__(args)
        self.filter_ma: list[str] = args.filter_monoamines

        self.src_col = 0
        self.tar_col = 1
        self.ma_col = 2
        self.receptor_col = 3

    def load(self, *args):
        csv_path = args[0]
        df = pd.read_csv(csv_path, header=None)

        self.logger.info(f'Filtering Neurons with monoamines: {self.filter_ma}')
        df = df[df[self.ma_col].isin(self.filter_ma)]

        if not len(df):
            raise Exception('Filtering results with an empty data')

        src_neurons_names = df.iloc[:, self.src_col]
        tar_neurons_names = df.iloc[:, self.tar_col]

        self.network.neuron_names = list(set(src_neurons_names) | set(tar_neurons_names))
        # make it deterministic between runs
        self.network.neuron_names.sort()
        neurons_indices = {ss: i for i, ss in enumerate(self.network.neuron_names)}

        for idx, (v1, v2, ma, receptor) in df.iterrows():
            if not self.args.allow_self_loops and v1 == v2:
                continue

            self.network.participating_neurons.add(v1)
            self.network.participating_neurons.add(v2)

            n1 = neurons_indices[v1]
            n2 = neurons_indices[v2]
            # TODO: some edges twice due to different receptors - increase the weight?
            # if self.network.graph.has_edge(n1, n2):
            #     print(v1, v2)
            self.network.graph.add_edge(n1, n2, monoamine=ma, receptor=receptor, polarity=None)

        self.network.use_monoamines = True


