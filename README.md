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
- **Triadic Census** - (networkx implementation) Vladimir Batagelj and Andrej Mrvar, "A subquadratic triad census algorithm for large sparse networks with small maximum degree" 2001
- **Direct** graph counting of specific patterns

##
Implemented random networks algorithms:
 - **Simple Markov-Chain**: R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).
 - **Erdos Renyi**: with probability p such that the average |E| of all random networks ~= |E| of the original network.  [(Wikipedia)](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model)
 - **Barabási–Albert model**: with m (# of edges to attach) = |E| / |N|. then random direction per edge is chosen. [(Wikipedia)](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model)
   
##
Supporting the following network formats:
- **simple** txt format: (v1, v2, w) per line
- **worm wiring** xlsx format: https://wormwiring.org/pages/adjacency.html
- **polarity** xlsx format from the paper: Fenyves BG, Szilágyi GS, Vassy Z, Sőti C, Csermely P. "Synaptic polarity and sign-balance prediction using gene expression data in the Caenorhabditis elegans chemical synapse neuronal connectome network"
- **Durbin** 1986 - https://www.wormatlas.org/neuronalwiring.html
