## Prokariotic Clustering based on Genome and Sequence Identity/Similarity Matrices

### ProKlust (“Prokaryotic Clusters”) was written with a focus on taxonomical data. It obtains, filters and visualizes clusters from multiple identity/similarity matrices using maximal clique enumeration (MCE) with settable cut-off points.

- The input pairwise identity/similarity matrix is formatted into triangular matrix using the average. Then, the matrix is formatted using the cut-off values chosen by the user, replacing values according to the criteria. If more than one matrix is used as input, each of the generated logical matrices (boolean) is multiplied. Afterwards, the igraph object is formed by the connection the nodes (for example, genomes or genes) which present the value 1. It uses a modified Bron–Kerbosch algorithm on “igraph”. 
- ProKlust allows the modification of graph data with the use of specific filters (or a combination/intersection of them) and renaming of nodes. It generates relevant cluster data as output and allows the possibility of graph visualization using the networkD3 library.

### How to install

- Install the devtools package: install.packages("devtools")
- Load devtools: library(devtools)
- Install this package: install_github("camilagazolla/ProKlust")
- Done! You can then load it by: library(ProKlust).
