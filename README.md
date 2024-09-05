Analyzing Network Motifs with Python.
<br/>
🪱 We focus on **C.elegans** connectome networks but the code can be applied to any network. 🪱

##
Network Motif discovery main steps:
1. Extract all possible subgraphs of a given size k
2. Enumerate the frequency of each subgraph
3. Test for a statistical significance of each subgraph using random networks.

##
Implemented enumeration algorithms:
- **Mfinder**, both an induced and non-induced version. R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002"
- **FANMOD** - S. Wernicke, “Efficient detection of network motifs” 2006
- **Triadic Census** - (a wrapper of networkx implementation) Vladimir Batagelj and Andrej Mrvar, "A subquadratic triad census algorithm for large sparse networks with small maximum degree" 2001

##
Implemented Large Motif search algorithms:
- **SIM** - Single Input moudle. Shai S. Shen-Orr, Ron Milo, Shmoolik Mangan & Uri Alon, "Network motifs in the transcriptional regulation  network of Escherichia coli"

##
Implemented random networks algorithms:
 - **Simple Markov-Chain**: R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
 - **Erdos Renyi**: with probability p such that the average |E| of all random networks ~= |E| of the original network.  [(Wikipedia)](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model)
 - **Barabási–Albert model**: with m (# of edges to attach) = |E| / |N|. then random direction per edge is chosen. [(Wikipedia)](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model)
 - **Anatomical constrained Markov-Chain**: allow switches based on the C.elegans nerve-ring distances, based on Zaslaver et al., 2022: "The synaptic organization in the Caenorhabditis elegans neural network suggests significant local compartmentalized computations"
##
Supporting the following network formats:
- **simple** txt format: (v1, v2, w) per line
- **worm wiring** format: https://wormwiring.org/pages/adjacency.html
- **polarity** format from the 2020 paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P. "Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical synapse neuronal connectome network"
- **Durbin** 1986 - https://www.wormatlas.org/neuronalwiring.html
- **multilayer connectome** format from the 2016 paper: Bentley B, Branicky R, Barnes CL, Chew YL, Yemini E, et al. "The Multilayer Connectome of Caenorhabditis elegans"

##
Post Motif-Search analysis' features:
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

