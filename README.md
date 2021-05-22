### ProKlust (“Prokaryotic Clusters”) was written with a focus on taxonomical data. It obtains, filters and visualizes clusters from multiple identity/similarity matrices.

ProKlust could be employed to analyze any identity/similarity matrix, such as ANI or barcoding gene identity. Additionally, it contains useful filter options to deal with taxonomical data. 

### How to install
```
library(devtools)
install_github("camilagazolla/ProKlust") # Install this package
library(ProKlust)
```

### Package dependencies
- [igraph](https://cran.r-project.org/web/packages/igraph/index.html)
- [networkD3](https://cran.r-project.org/web/packages/networkD3/index.html)

### How to use

- Use function prokluster to obtain the clusters.
- Use function plotc to plot.

#### Inputs
- Obligatory tabbed-delimited pairwise identity/similarity matrix(ces), such as the ones generated with pyANI.

**IMPORTANT:** if the user wishes to send a list of matrices (instead of a vector of file names), he **MUST** convert it to a list, e.g.
```
percentage <- read.table(file = "ANIb_percentage_identity.tab", header = T, row.names = 1, sep = "\t")
coverage <- read.table(file = "ANIb_alignment_coverage.tab", header = T, row.names = 1, sep = "\t")
filesList <- list(percentage, coverage)
thresholds <- c(0.95, 0.70)
basicResult <- prokluster(files = filesList, cutoffs = thresholds)
plotc(basicResult$graph)
```
- Optional annotation table text file with a specific format: each line must be "(previous name)<Tab>(new name)<New Line>", where (previous name) and (new name) can be alphanumeric and some special characters.

#### Filters
- filterRemoveIsolated: remove isolated nodes i.g. nodes that do not form groups/clusters; 
- filterRemoveLargerComponent and filterOnlyLargerComponent: remove or retain only the component containing the highest number of nodes;
- filterDifferentNamesConnected: retain groups of connected nodes containing more than one binomial species name;
- filterSameNamesNotConnected: retain groups of unconnected nodes containing the same species names.

#### Outputs

- maxCliques: the maximal clique is the largest subset of nodes in which each node is directly connected to every other node in the subset; 
- components: contains the isolated nodes or groups formed of complete graphs; 
- graph: an igraph object graph, that can be further handled by the user;
- plot: where the final graph could be promptly visualized with forceNetwork function from the networkD3 R package.

A genome/gene that is part of a component does not necessarily share identity/similarity values above the established cut-off with all the other genomes/genes of that component, but it must share an identity/similarity value above the cut-off for at least one other genome/gene. Cliques, instead, are formed by genome/gene that all share identity/similarity values above the chosen criteria. A genome/gene could belong at the same time to different cliques within the same component.

#### Examples
```
#Example 1.1
basicResult1.1 <- prokluster(files = "ANIb_percentage_identity.tab", cutoffs = 0.9)
basicResult1.1
plotc(basicResult1.1$graph)

#Example 1.2
percentage <- read.table(file = "ANIb_percentage_identity.tab", header = T, row.names = 1, sep = "\t")
basicResult1.2 <- prokluster(files = percentage, cutoffs = 0.9)

#Example 2.1
files <- c("ANIb_percentage_identity.tab", "ANIb_alignment_coverage.tab")
thresholds <- c(0.95, 0.70)
renamedResults1.1 <- prokluster(files = files, cutoffs = thresholds, nodesDictionary = "dictionary.tab", filterRemoveIsolated = TRUE)

#Example 2.2
coverage <- read.table(file = "ANIb_alignment_coverage.tab", header = T, row.names = 1, sep = "\t")
filesList <- list(percentage, coverage)
basicResult2.2 <- prokluster(files = filesList, cutoffs = thresholds)

#Example 3
renamedResults2 <- prokluster(files = files, cutoffs = thresholds, nodesDictionary = "dictionary.tab", filterDifferentNamesConnected = TRUE)

#Example 4
nodesNames <- read.table(file= "dictionary.tab", sep = "\t", header = F, stringsAsFactors=FALSE)
renamedResults3 <- prokluster(files = files, cutoffs = thresholds, nodesPreviousNames = nodesNames$V1, nodesTranslatedNames = nodesNames$V2, filterSameNamesNotConnected = T)
```

#### Example for [FastANI](https://github.com/ParBLiSS/FastANI):
  
```
$ cd bins #dir with genomes
$ mkdir out
$ ls *fna > list
$ mv list out
$ for f in *fna; do fastANI -q "${f}" --rl out/list -o "${f}.fastANI" --minFraction 0; mv "${f}.fastANI" out; done
$ cd out/
$ cat *ANI > fastANIout.txt
```

Generating a tabbed-delimited "pairwise" identity matrix on R:
```
library(ProKlust)
library(tidyr)

# Importing fastANI results
identity <- (read.table(file = "fastANIout.txt", sep = "\t")) [1:3]
identity <- pivot_wider(identity, names_from =V1, values_from = V3)
identity <- as.data.frame(identity)
rownames(identity) <- identity$V2
identity.sorted <- identity[order(identity["V2"]),]
identity.sorted[,1] <- NULL

basicResult <- prokluster(file = identity.sorted, cutoffs = 95)
basicResult
plotc(basicResult$graph)
```
### Workflow


<img src="https://user-images.githubusercontent.com/64544051/109963334-0d672180-7ccb-11eb-85f5-238f9cefbe76.JPEG" width="50%" height="50%">


A) The average of each pair from the pairwise input matrix/matrices is/are obtained. A Boolean matrix/matrices is/are obtained according to the cut-off values chosen by the user. If more than one matrix is used as input, the final generated matrix is obtained by multiplying the elements of the matrices. A graph is formed by connecting the nodes which present the positive values. In this example, nodes correspond to genomes and edges correspond to ANI ≥ 95% with coverage alignment ≥ 50%. The data could be filtered to retain components containing more than one species name or unconnected nodes containing the same species names. 

B) Overview of the hierarchical-based clustering approach. These approaches return tree-shaped diagrams with non-overlapping clusters.

### Citation
If you use ProKlust in your research please cite:

>Volpiano CG, Sant’Anna FH, Ambrosini A, de São José JFB, Beneduzi A, Whitman WB, de Souza EM, Lisboa BB, Vargas LK and Passaglia LMP (2021) Genomic Metrics Applied to _Rhizobiales_ (_Hyphomicrobiales_): Species Reclassification, Identification of Unauthentic Genomes and False Type Strains. Front. Microbiol. 12:614957. https://doi.org/10.3389/fmicb.2021.614957
