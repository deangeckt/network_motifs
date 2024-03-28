import numpy as np
import pandas as pd


def common_adj_mat_formatting(adj_mat: np.ndarray, neuron_names: list[str], output_name: str):
    participating_neurons = set()
    adj_output_file_name = f'{output_name}_adj.txt'
    with open(adj_output_file_name, "w") as f:
        amount_of_synapses = 0
        amount_of_edges = 0
        for i in range(len(adj_mat)):
            for j in range(len(adj_mat)):
                synapses = adj_mat[i, j]
                if not np.isnan(synapses):
                    amount_of_synapses += synapses
                    amount_of_edges += 1
                    f.write(f'{i} {j} {synapses}\n')
                    participating_neurons.add(i)
                    participating_neurons.add(j)

    neuron_output_file_name = f'{output_name}_neurons.txt'
    with open(neuron_output_file_name, "w") as f:
        [f.write(f'{neuron}\n') for neuron in neuron_names]

    print(f'total number of neurons: {len(neuron_names)}')
    print(f'total number of synapses: {amount_of_synapses}')
    print(f'total number of network edges: {amount_of_edges}')
    print(f'total number of network nodes: {len(participating_neurons)}')

    print(f'saved adj matrix file to: {adj_output_file_name}')
    print(f'saved neurons names file to: {neuron_output_file_name}')


def convert_worm_wiring_xlsx(xlsx_path: str, sheet_name: str, output_name: str, row_start: int, col_start: int,
                             amount: int):
    """
    https://wormwiring.org/pages/adjacency.html xlsx format
    convert to the simple txt format: (v1, v2, w) per line.
    convert to the simple txt format: (v1, v2, w) per line. w isn't the weight but the number of synapses
    """
    xls = pd.ExcelFile(xlsx_path)
    df = xls.parse(sheet_name, header=None)

    names_idx = 2

    adj_mat = df.iloc[row_start:row_start + amount, col_start:col_start + amount]
    neurons_col = list(df.iloc[row_start:row_start + amount, names_idx])
    neurons_row = list(df.iloc[names_idx, col_start:col_start + amount])
    if not neurons_row == neurons_col:
        print(f'invalid xlsx file or indices')
        return

    adj_mat = np.array(adj_mat)
    common_adj_mat_formatting(adj_mat=adj_mat, neuron_names=neurons_col, output_name=output_name)


def convert_durbin_txt(file_path: str, output_name: str, filter_syn_type: str, filter_recon: str, filter_joint: bool):
    """
    https://www.wormatlas.org/neuronalwiring.html - Neuronal Connectivity I: by R. Durbin
    convert to the simple txt format: (v1, v2, w) per line. w isn't the weight but the number of synapses
    : param filter_syn_type: either 'chem', 'gap'
    : param filter_recon: either: one of 'JSH', 'N2U'
    : param filter_joint: relevant for 'chem' synapse only - either True or False
    """
    neurons_names = set()
    data = []

    with (open(file_path, "r") as f):
        for line in f.readlines():
            line = tuple(line.strip().split())
            if len(line) == 5:
                n1, n2, synapse_type, reconstruction, num_of_synapses = line
            elif len(line) == 4:
                n1, n2, synapse_type, num_of_synapses = line
                reconstruction = 'none'
            else:
                print('invalid file')
                return
            data.append((n1, n2, synapse_type, reconstruction, int(num_of_synapses)))

    data = list(filter(lambda x: x[3] == filter_recon, data))
    gap_data = list(filter(lambda x: x[2] == 'Gap_junction', data))
    chem_data = list(filter(lambda x: x[2] != 'Gap_junction', data))

    if filter_joint:
        chem_data = list(filter(lambda x: 'joint' in x[2], chem_data))
    else:
        chem_data = list(filter(lambda x: 'joint' not in x[2], chem_data))

    if filter_syn_type == 'gap':
        data = gap_data
    elif filter_syn_type == 'chem':
        data = chem_data
    else:
        print(f'invalid syn type: {filter_syn_type}')
        return

    for n1, n2, _, _, _ in data:
        neurons_names.add(n1)
        neurons_names.add(n2)

    neurons = {n: i for i, n in enumerate(list(neurons_names))}
    N = len(neurons)
    adj_mat = np.empty((N, N))
    adj_mat.fill(np.nan)

    for n1, n2, synapse_type, _, num_of_synapses in data:
        v1 = neurons[n1]
        v2 = neurons[n2]

        if 'Receive' in synapse_type:
            v1 = neurons[n2]
            v2 = neurons[n1]
        if not np.isnan(adj_mat[v1, v2]):
            adj_mat[v1, v2] += num_of_synapses
        else:
            adj_mat[v1, v2] = num_of_synapses

    common_adj_mat_formatting(adj_mat=adj_mat, neuron_names=list(neurons_names), output_name=output_name)


if __name__ == "__main__":
    # convert_worm_wiring_xlsx('../networks/Cook_2019/SI 2 Synapse adjacency matrices.xlsx',
    #                          'herm chem synapse adjacency',
    #                          '2020_si_2_herm_chem_synapse', 60, 60, 272)
    #
    # convert_worm_wiring_xlsx('../networks/Cook_2019/SI 2 Synapse adjacency matrices.xlsx',
    #                          'herm gap jn synapse adjacency',
    #                          '2020_si_2_herm_gap', 60, 60, 272)
    #
    # convert_worm_wiring_xlsx('../networks/Cook_2019/SI 5 Connectome adjacency matrices, corrected July 2020.xlsx',
    #                          'hermaphrodite chemical',
    #                          '2020_si_5_herm_chem_synapse', 23, 53, 272)
    convert_durbin_txt("../networks/Durbin_1986/neurodata.txt",
                       'durbin_herm_chem_synapses',
                       'chem',
                       'N2U',
                       False)
