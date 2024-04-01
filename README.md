Analyzing Network Motifs with Python.
<br/>
🪱 We focus on **C.elegans** networks but the code can be applied to any network. 🪱

##
Network Motif discovery main steps:
1. Extract all possible subgraphs of a given size k
2. Enumerate the frequency of each subgraph
3. Test for a statistical significance of each subgraph using random networks.

##
Implemented enumeration algorithms:
- **Mfinder**, both an induced and non-induced version ( R. Milo, S. Shen-Orr, S. Itzkovitz, N. Kashtan,D. Chklovskii, and U. Alon, “Network motifs: simple building blocks of complex networks.” Science, vol. 298, no. 5594, pp. 824–827, October 2002.)
- **Direct** graph counting of specific patterns

##
Implemented random networks algorithms:
 - **Simple Markov-Chain**: R. Kannan, P. Tetali, S. Vempala, Random Struct. Algorithms 14, 293 (1999).

##
Network format:
<br/>
We work with a simple network format - a txt file where each line is a tupple: (v1, v2, w), spereted by spaces.
