
# Detect and analyze Network Motifs of any network, with Python. #
🪱 We focus on **_C.elegans_** connectome networks but the code can be applied to any network. 🪱

## We implemented and support: ##

#### Enumeration algorithms: ####
- **Mfinder**, both an induced and non-induced version. R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002"
- **FANMOD** - S. Wernicke, “Efficient detection of network motifs” 2006
- **Triadic Census** - (a wrapper of networkx implementation) Vladimir Batagelj and Andrej Mrvar, "A subquadratic triad census algorithm for large sparse networks with small maximum degree" 2001

#### Large Motif search algorithms: ####
- **SIM** - Single Input moudle. Shai S. Shen-Orr, Ron Milo, Shmoolik Mangan & Uri Alon, "Network motifs in the transcriptional regulation  network of Escherichia coli"

#### Random networks algorithms: ####
 - **Simple Markov-Chain**: R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
 - **Erdos Renyi**: with probability p such that the average |E| of all random networks ~= |E| of the original network.  [(Wikipedia)](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model)
 - **Barabási–Albert model**: with m (# of edges to attach) = |E| / |N|. then random direction per edge is chosen. [(Wikipedia)](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model)
 - **Anatomical constrained Markov-Chain**: allow switches based on the _C.elegans_ nerve-ring distances, based on Zaslaver et al., 2022: "The synaptic organization in the Caenorhabditis elegans neural network suggests significant local compartmentalized computations"

#### Network formats: ####
- **simple** txt format: (v1, v2, w) per line
- **worm wiring** format: https://wormwiring.org/pages/adjacency.html
- **polarity** format from the 2020 paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P. "Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical synapse neuronal connectome network"
- **Durbin** 1986 - https://www.wormatlas.org/neuronalwiring.html
- **multilayer connectome** format from the 2016 paper: Bentley B, Branicky R, Barnes CL, Chew YL, Yemini E, et al. "The Multilayer Connectome of Caenorhabditis elegans"

## Usage: ##

This repository contains two main parts: first is a generic network-motif detection software, second is a post-motif-search analysis features spread out in jupyter notebooks.

#### Motif detection software: ####
The motif_search_main.py is the main python script which can be ran via any shell / script in a command-line-argument fashion (it also includes the documentations):

Example from the terminal:
```
motif_search_main.py --help

optional arguments:
  -h, --help            show this help message and exit
  -rs RANDOM_SEED, --random_seed RANDOM_SEED
                        random seed for the entire program
  -lf LOG_FILE, --log_file LOG_FILE
                        file path to save log results
  -bf BIN_FILE, --bin_file BIN_FILE
                        file path to save binary results
  -it {simple_adj_txt,worm_wiring_xlsx,polarity_xlsx,durbin_txt,multilayer,graph,binary_network_file}, --input_type {simple_adj_txt,worm_wiring_xlsx,polarity_xlsx,durbin_txt,multilayer,graph,binary_network_file}
                        the type of the input network
  -inf INPUT_NETWORK_FILE, --input_network_file INPUT_NETWORK_FILE
                        file path of the input network
  -ing INPUT_NETWORK_GRAPH [INPUT_NETWORK_GRAPH ...], --input_network_graph INPUT_NETWORK_GRAPH [INPUT_NETWORK_GRAPH ...]
                        a graph: list of strings (tuples) where each is an edge. in the format: ["1 2" "2 3" ...]
  -rmc, --run_motif_criteria
                        run full motif search with motif criteria tests
  -sa {mfinder_i,mfinder_ni,fanmod,triadic_census,specific}, --sub_graph_algorithm {mfinder_i,mfinder_ni,fanmod,triadic_census,specific}
                        sub-graph enumeration algorithm
  -k K, --k K           the size of sub-graph / motif to search in the enumeration algorithm
  -sim SIM, --sim SIM   the maximum size of control size in the SIM search algorithm
  -uim, --use_isomorphic_mapping
                        run (pre motif search) isomorphic sub-graphs search
  -asl, --allow_self_loops
                        allow self loops in the (pre motif search) isomorphic sub-graphs search
  -st SYNAPSE_THRESHOLD, --synapse_threshold SYNAPSE_THRESHOLD
                        filter neurons with >= # synapses (only in neuron networks files)
  -fsy {chem,gap,all}, --filter_syn_type {chem,gap,all}
                        filter synapse type, supported in durbin and worm_wiring networks
  -fsx {herm,male}, --filter_sex_type {herm,male}
                        filter sex type, supported in durbin and worm_wiring networks
  -fnrn, --filter_nerve_ring_neurons
                        filter the neuronal connectome with only the (180) nerve ring neurons
  -fp {+,-,no pred,complex} [{+,-,no pred,complex} ...], --filter_polarity {+,-,no pred,complex} [{+,-,no pred,complex} ...]
                        polarity: filter neurons with polarity
  -fpn {Glu,GABA,ACh,0} [{Glu,GABA,ACh,0} ...], --filter_prim_nt {Glu,GABA,ACh,0} [{Glu,GABA,ACh,0} ...]
                        polarity: filter neurons with primary neurotransmitter
  -fma {dopamine,octopamine,serotonin,tyramine} [{dopamine,octopamine,serotonin,tyramine} ...], --filter_monoamines {dopamine,octopamine,serotonin,tyramine} [{dopamine,octopamine,serotonin,tyramine} ...]
                        Monoamines: filter neurons with MA transmitter
  -r {markov_chain,nerve_ring_markov_chain,erdos_renyi,barabasi}, --randomizer {markov_chain,nerve_ring_markov_chain,erdos_renyi,barabasi}
                        main randomizer algorithm in a full motif search
  -na NETWORK_AMOUNT, --network_amount NETWORK_AMOUNT
                        amount of random networks to generate in a full motif search
  -sf SWITCH_FACTOR, --switch_factor SWITCH_FACTOR
                        number of switch factors done by the markov chain randomizer
  -a ALPHA, --alpha ALPHA
                        motif criteria alpha for testing p value significance
  -ft FREQUENCY_THRESHOLD, --frequency_threshold FREQUENCY_THRESHOLD
                        motif criteria frequency threshold test
  -ut UNIQUENESS_THRESHOLD, --uniqueness_threshold UNIQUENESS_THRESHOLD
                        motif criteria uniqueness threshold test
  -uut, --use_uniq_criteria
                        whether to use the uniqueness test
```

Example via .bat script:
```
@echo off
python motif_search_main.py ^
--input_type durbin_txt ^
--input_network_file "networks/data/Durbin_1986/neurodata.txt" ^
--synapse_threshold 10 ^
--filter_sex_type herm ^
--filter_syn_type chem ^
--sub_graph_algorithm mfinder_i ^
--k 3 ^
--sim 1 ^
--bin_file "durbin_herm_chem_k4_m5.bin" ^
--run_motif_criteria
```

In both cases it is mandatory to run from the **root** folder.

Basic output example:
```
Filtering Neurons of: herm
Filtering Synapses of type: chem
Network file name: SI 2 Synapse adjacency matrices.xlsx

Network properties:
        Neurons in the network: 272
        Participating Neurons (in the graph): 250

        Participating Nodes are neurons with at least: 15 synapses
        Synapses in the network: 9580
        Synapses in the graph: 1587
        Gaps in the network: 0
        Gaps in the graph: 0
        Nodes: 87
        Edges: 80
        Average clustering coefficient: 0.024
        Average shortest path (undirected): 1.415
        Density: 0.011
        Degree: Mean: 1.839 Std: 1.66 Median: 1.0 Max: 12 (node: RIAR)
        In-Degree: Mean: 0.92 Std: 1.147 Median: 1.0 Max: 5 (node: AVAR)
        Out-Degree: Mean: 0.92 Std: 1.127 Median: 1.0 Max: 8 (node: RIAR)

Motif criteria:
        alpha: 0.01
        use uniqueness: False
        uniqueness threshold: 3
        frequency threshold: 0.1

Sub graph enumeration algorithm: fanmod
Search using k: 3
Using isomorphic mapping: True
SIM sub graph search using max-control size: 1
Allow self loops: False
Using nerve ring neurons only: False
Full motif search: True
Sub Graph search timer [Sec]: 0.07
SIM search timer [Sec]: 0.0
Populate polarity motifs timer [Sec]: 0.0

Randomizer: using markov_chain algorithm
Randomizer: generating 1000 random networks
Markov chain switch factor: 10
Markov chain iterations: 800
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1000/1000 [00:11<00:00, 86.29it/s] 
Markov chain success ratio: 0.90157
average # edges of all random networks: 80.0
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1000/1000 [01:05<00:00, 15.22it/s] 

Motif candidates found: 7
Total number of sub graphs found: 154
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  ID  |     Name     |  Adj_mat  |  is-motif  |  N_real  |  N_rand  |  Std  |  Z_score  |  uniq  |  is-significant  |  is-freq  |  is-uniq  |  is-anti-freq  |
+======+==============+===========+============+==========+==========+=======+===========+========+==================+===========+===========+================+
|  6   |   fan_out    | [[0 1 1]  | anti-motif |    35    |  47.73   | 4.78  |   -2.66   |   3    |       True       |   False   |    n/a    |      True      |
|      |              |  [0 0 0]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  12  |   cascade    | [[0 0 1]  | anti-motif |    52    |  73.39   | 8.51  |   -2.51   |   1    |       True       |   False   |    n/a    |      True      |
|      |              |  [1 0 0]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  14  |     n/a      | [[0 1 1]  |    none    |    12    |   2.91   | 4.45  |   2.04    |   1    |      False       |   True    |    n/a    |     False      |
|      |              |  [1 0 0]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  36  |    fan_in    | [[0 0 1]  | anti-motif |    44    |  50.91   | 2.91  |   -2.37   |   8    |       True       |   False   |    n/a    |      True      |
|      |              |  [0 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  38  | feed_forward | [[0 1 1]  |   motif    |    5     |   1.15   | 1.07  |   3.59    |   1    |       True       |   True    |    n/a    |     False      |
|      |              |  [0 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  46  |     n/a      | [[0 1 1]  |    none    |    0     |   0.06   | 0.29  |   -0.21   |   0    |      False       |   False   |    n/a    |      True      |
|      |              |  [1 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [0 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  74  |     n/a      | [[0 1 0]  |    none    |    5     |   1.76   | 2.62  |   1.23    |   1    |      False       |   True    |    n/a    |     False      |
|      |              |  [1 0 0]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  78  |     n/a      | [[0 1 1]  |    none    |    0     |   0.04   | 0.25  |   -0.15   |   0    |      False       |   False   |    n/a    |      True      |
|      |              |  [1 0 0]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
|  98  |     n/a      | [[0 1 0]  |    none    |    1     |   0.19   | 0.42  |   1.92    |   1    |      False       |   True    |    n/a    |     False      |
|      |              |  [0 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
| 102  |     n/a      | [[0 1 1]  |    none    |    0     |   0.03   | 0.16  |   -0.17   |   0    |      False       |   False   |    n/a    |      True      |
|      |              |  [0 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
| 108  |     n/a      | [[0 0 1]  |    none    |    0     |   0.02   | 0.14  |   -0.14   |   0    |      False       |   False   |    n/a    |      True      |
|      |              |  [1 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
| 110  |     n/a      | [[0 1 1]  |    none    |    0     |   0.00   | 0.04  |   -0.04   |   0    |      False       |   False   |    n/a    |      True      |
|      |              |  [1 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 0 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
| 238  |     n/a      | [[0 1 1]  |    none    |    0     |   0.00   | 0.00  |   0.00    |   0    |      False       |   False   |    n/a    |     False      |
|      |              |  [1 0 1]  |            |          |          |       |           |        |                  |           |           |                |
|      |              |  [1 1 0]] |            |          |          |       |           |        |                  |           |           |                |
+------+--------------+-----------+------------+----------+----------+-------+-----------+--------+------------------+-----------+-----------+----------------+
```


#### Post Motif-Search analysis' features: ####

These utilities are located under /post-motif-analysis folder, they need a binary output file (from the previous step) to work:
```
--bin_file "your_results_file.bin"
```

- Plot Motifs's Z-score vs N_real:
 ![image](https://github.com/user-attachments/assets/c29c8594-4330-42cb-b6e0-a2b9c1630f1e)

- Compare multiple networks normalized and discretized Z-scores:
  ![image](https://github.com/user-attachments/assets/203ab81b-2843-48d8-a1d7-476ea93a88d2)
![image](https://github.com/user-attachments/assets/d92747bf-21fa-4428-b267-8e3b18ffe8dd)
![image](https://github.com/user-attachments/assets/2c47ea4d-de5d-4a48-b46e-2c376f631775)


- Compare 2 networks Z-scores (via a scatter):
![image](https://github.com/user-attachments/assets/6e2c6c63-7163-4684-8e2f-1b34b29c88c0)

- Sort Motifs and Anti-Motifs based on Z-score:
  
 ![image](https://github.com/deangeckt/network_motifs/assets/24900065/730c2962-0045-4a82-beea-bc1b86ecbc06)
 ![image](https://github.com/deangeckt/network_motifs/assets/24900065/35ed62ef-d0a1-4c4a-9f69-00aa88a9bde3)

- Analyze the frequency of nodes per their role in a given Motif:
![image](https://github.com/deangeckt/network_motifs/assets/24900065/c4092a87-0bbe-4994-b2e1-13a03312a979)

- Analyze the frequency of specific nodes in all Motifs:
![image](https://github.com/deangeckt/network_motifs/assets/24900065/a62d074c-a5fd-46ec-810a-1d999ab714b8)

- Filter (and later draw) all **subgraphs** where node N was part of motif M in role R, e.g.,:
![image](https://github.com/user-attachments/assets/601ceb0a-3644-4950-a57e-6f8f044ac912)

- Rank any set of subgraphs by different metrices (degree, clustering coefficient)
![image](https://github.com/user-attachments/assets/489a80ed-e52e-4a84-9a03-51075568e119)

- Explore larger-than-k=3 (beasts) Moitfs:
  ![image](https://github.com/user-attachments/assets/14afc5c1-d425-4a57-8374-7e1d923375d2)

