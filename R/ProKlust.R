# Make matrix in 'boolean' format
makeMatrixBoolean <- function(matrix, cutoff)
{
  matrix_upper <- matrix
  matrix_upper[upper.tri(matrix_upper)] <- NA
  matrix_lower <- matrix
  matrix_lower[lower.tri(matrix_lower)] <- NA

  # medianMatrix holds the mean value of upper and lower
  matrix_lower <- t(matrix_lower)
  medianMatrix <- ((matrix_upper+ matrix_lower)/2)

  # replace medianMatrix to boolean values according to cutoff
  medianMatrix[medianMatrix < cutoff] <- 0 #replace values lower than cutoff with 0
  medianMatrix[medianMatrix >= cutoff] <- 1 #replace values higher or equal to cutoff with 1

  return(medianMatrix)
}

# Reflect the matrix
mirrorMatrix <- function(matrix)
{
  matrixT <- t(matrix)
  mirroredMatrix <- matrix
  mirroredMatrix[,] <- NA #same size
  mirroredMatrix[upper.tri(mirroredMatrix)] <- matrixT[upper.tri(matrixT)]
  mirroredMatrix[lower.tri(mirroredMatrix)] <- matrix[lower.tri(matrix)]
  return(mirroredMatrix)
}

# Open matrix files and compute cutoffs to combine them
matricesCombination <- function(filenames, cutoffs)
{
  if (length(filenames) == length(cutoffs) && !is.na(filenames) && !is.na(cutoffs))
  {
    if (is.na(length(filenames)))
    {
      currentMatrix <- read.table(file = filenames, header = T, row.names = 1, sep = "\t")
      currentBooleanMatrix <- makeMatrixBoolean(currentMatrix, cutoffs)
      return(mirrorMatrix(currentBooleanMatrix))
    }

    for (i in 1:length(filenames))
    {
      currentMatrix <- read.table(file = filenames[i], header = T, row.names = 1, sep = "\t")
      currentBooleanMatrix <- makeMatrixBoolean(currentMatrix, cutoffs[i])

      if (i == 1)
      {
        matrixMultiplication <- currentBooleanMatrix
      }
      else
      {
        matrixMultiplication <- matrixMultiplication * currentBooleanMatrix
      }
    }

    return(mirrorMatrix(matrixMultiplication))
  }
  else
  {
    stop("filenames and cutoffs must have the same number of elements.")
  }
}

# Filter igraph to remove groups containing only one species name
filterDifferentNamesConnected <- function(g)
{
  repetition <- TRUE

  while(repetition == TRUE)
  {
    # get current df of g
    gcomponents <- as.data.frame(components(g)$membership)
    gcomponents$spp <-rownames(gcomponents)
    gcomponents$spp <- gsub("(\\w+\\s\\w+).*", "\\1", gcomponents$spp)
    gcomponents <- gcomponents[order(gcomponents$`components(g)$membership`),]

    # generate df with repeated groups
    repeatedGroups <- gcomponents
    repeatedNames <- table(repeatedGroups$`components(g)$membership`) > 1
    repeatedGroups <- subset(repeatedGroups, `components(g)$membership` %in% names(repeatedNames)[repeatedNames])

    groupsToStay = vector()
    uniqueGroups <- unique(repeatedGroups[,1])

    for (uniqueGroup in uniqueGroups)
    {
      repeatedData <- repeatedGroups[repeatedGroups$`components(g)$membership` == uniqueGroup,]
      uniqueNames <- unique(repeatedData[,2])

      # if there is a name with more than one group, it should stay
      if (length(uniqueNames) > 1)
      {
        groupsToStay <- c(groupsToStay, uniqueGroup)
      }
    }

    groupsToRemove <- gcomponents$`components(g)$membership`[!(gcomponents$`components(g)$membership` %in% groupsToStay)]
    groupsToRemove <-  unique(groupsToRemove)

    if (length(groupsToRemove) == 1)
    {
      repetition <-  FALSE
    }

    g <- induced_subgraph(g, V(g)[components(g)$membership != groupsToRemove[1]])
  }

  return(g)
}

# Filter igraph groups that contain same name and are not connected.
filterSameNamesNotConnected <- function(g)
{
  repetition <- TRUE

  while(repetition == TRUE)
  {
    # get current df of g
    gcomponents <- as.data.frame(components(g)$membership)
    gcomponents$spp <-rownames(gcomponents)
    gcomponents$spp <- gsub("(\\w+\\s\\w+).*", "\\1", gcomponents$spp)
    gcomponents <- gcomponents[order(gcomponents$`components(g)$membership`),]

    # generate df with repeated spp
    repeatedNames <- gcomponents
    repeatedGroups <- table(repeatedNames$spp) > 1
    repeatedNames <- subset(repeatedNames, spp %in% names(repeatedGroups)[repeatedGroups])

    groupsToStay = vector()
    uniqueNames <- unique(repeatedNames[,2])

    for (uniqueName in uniqueNames)
    {
      repeatedData <- repeatedNames[repeatedNames$spp == uniqueName,]
      uniqueGroups <- unique(repeatedData[,1])

      # if there is a name with more than one group, it should stay
      if (length(uniqueGroups) > 1)
      {
        groupsToStay <- c(groupsToStay, uniqueGroups)
      }
    }

    groupsToRemove <- gcomponents$`components(g)$membership`[!(gcomponents$`components(g)$membership` %in% groupsToStay)]
    groupsToRemove <-  unique(groupsToRemove)

    if (length(groupsToRemove) == 1)
    {
      repetition <-  FALSE
    }

    g <- induced_subgraph(g, V(g)[components(g)$membership != groupsToRemove[1]])
  }

  return(g)
}

# Process filters on the original graph. If more than one filter is chosen, the resulted graph will try to be an intersection of the filters.
runGraphFilters <- function(g, filterRemoveIsolated, filterRemoveLargerComponent, filterOnlyLargerComponent, filterDifferentNamesConnected, filterSameNamesNotConnected)
{
  gFiltered <- g
  hasOneFilter <-  FALSE

  if (filterRemoveIsolated)
  {
    Isolated = (degree(g)== 0)
    gFiltered <- delete.vertices(g, Isolated)
    hasOneFilter <-  TRUE
  }

  if (filterRemoveLargerComponent)
  {
    gfilterRemoveLargerComponent <- induced_subgraph(g, V(g)[components(g)$membership != which.max(components(g)$csize)])

    if (hasOneFilter)
    {
      gFiltered <- graph.intersection(gFiltered, gfilterRemoveLargerComponent, byname = "auto", keep.all.vertices = FALSE)
    }
    else
    {
      hasOneFilter <-  TRUE
      gFiltered <- gfilterRemoveLargerComponent
    }
  }

  if (filterOnlyLargerComponent)
  {
    gfilterOnlyLargerComponent <- induced_subgraph(g, V(g)[components(g)$membership == which.max(components(g)$csize)])

    if (hasOneFilter)
    {
      gFiltered <- graph.intersection(gFiltered, gfilterOnlyLargerComponent, byname = "auto", keep.all.vertices = FALSE)
    }
    else
    {
      hasOneFilter <-  TRUE
      gFiltered <- gfilterOnlyLargerComponent
    }
  }

  if (filterDifferentNamesConnected)
  {
    gfilterDifferentNamesConnected <- filterDifferentNamesConnected(g)

    if (hasOneFilter)
    {
      gFiltered <- graph.intersection(gFiltered, gfilterDifferentNamesConnected, byname = "auto", keep.all.vertices = FALSE)
    }
    else
    {
      hasOneFilter <-  TRUE
      gFiltered <- gfilterDifferentNamesConnected
    }
  }

  if (filterSameNamesNotConnected)
  {
    gfilterSameNamesNotConnected <- filterSameNamesNotConnected(g)

    if (hasOneFilter)
    {
      gFiltered <- graph.intersection(gFiltered, gfilterSameNamesNotConnected, byname = "auto", keep.all.vertices = FALSE)
    }
    else
    {
      hasOneFilter <-  TRUE
      gFiltered <- gfilterSameNamesNotConnected
    }
  }

  return(gFiltered)
}

#' Plot the computed clusters igraph.
#'
#' networkD3 library is used to plot the graph and a HTML copy file named network.html is created.
#'
#' @param igraph Obligatory igraph object.
#' @param plotFontSize Optional numeric font size in pixels for the node text labels.
#' @param plotLinkDistance Optional numeric or character string. Either numberic fixed distance between the links in pixels (actually arbitrary relative to the diagram's size). Or a JavaScript function, possibly to weight by Value. For example: linkDistance = JS("function(d){return d.value * 10}").
#' @param nodeAttraction Optional numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
#'
#' @return The actual plot of the network (an igraph which was transformed into graphical plot by networkD3).
#'
#' @examples
#' files <- c("ANIb_percentage_identity.tab", "ANIb_alignment_coverage.tab")
#' thresholds <- c(0.95, 0.70)
#' results <- prokluster(filenames = files, cutoffs = thresholds, nodesDictionary = "dictionary.tab")
#' plottedD3 <- plotc(results$graph)
#'
#' @export
plotc <- function(igraph, plotFontSize = 10, plotLinkDistance = 60, nodeAttraction = -10)
{
  # Transform Igraph format in something readable by networkD3
  library(networkD3)

  g_network=igraph_to_networkD3(igraph)

  if (nrow(g_network$links) == 0)
  {
    g_network$links[nrow(g_network$links) + 1, ] <- 0
  }

  g_network$nodes$group <-  1
  ColourScale <- 'd3.scaleOrdinal().domain(["1"]).range(["#320b51"]);'

  oldWarningOption <- getOption("warn")
  options(warn = -1)

  plottedNetwork <- forceNetwork(Links = g_network$links,
                                 Nodes = g_network$nodes,
                                 Source = 'source', Target = 'target',
                                 NodeID = 'name',
                                 Group ='group',
                                 opacity = 0.9,
                                 opacityNoHover = 1,
                                 colourScale = JS(ColourScale),
                                 zoom = T,
                                 linkWidth = 1,
                                 linkColour = "#e3dede",
                                 fontFamily = "arial",
                                 fontSize = plotFontSize, # size of the node names
                                 linkDistance = plotLinkDistance, # distance between node. Increase this value to have more space between nodes
                                 charge = nodeAttraction) # numeric value indicating either the strength of the node repulsion (negative value)

  options(warn = oldWarningOption)

  saveNetwork(plottedNetwork, "network.html")
  return(plottedNetwork)
}

#' Compute prokariotic clusters for the given matrices and cutoffs.
#'
#' Obtains and filters clusters of the identity/similarity matrices, identifying MCE by configuring settable cut-off points for each of the multiple matrices entries. Returns relevant graph information.
#'
#' @param filenames Obligatory tabbed-delimited pairwise identity/similarity matrix input file(s) name(s). Either a character vector or a vector of character vectors. If it's the second option, then the file name sequence must match the sequence of the cutoff list.
#' @param cutoffs Obligatory cutoff number for the given input file(s). Either a single number or a vector of numbers. If it's the second option, then the file name sequence must match the sequence of the cutoff list. Must be given special attention to format the cutoff with the given matrix.
#' @param nodesDictionary Optional annotation table text file with a specific format: each line must be "(previous name)<Tab>(new name)<New Line>", where (previous name) and (new name) can be alphanumeric and some special characters.
#' @param nodesPreviousNames Optional vector containing character vectors that represent the previous names (to be renamed) of the graph nodes. Must be of the same size as nodesTranslatedNames.
#' @param nodesTranslatedNames Optional vector containing character vectors that represent the new names of the graph nodes. Must be of the same size as nodesPreviousNames
#' @param filterRemoveIsolated Optional boolean parameter that allows the removal of isolated nodes (nodes that have no connection) of the graph. If more than one type of filter is chosen, the intersection of the filters (executed in the original graph) will be returned.
#' @param filterRemoveLargerComponent Optional boolean parameter that allows the removal of the larger (with the most connections) component/group of the graph. If more than one type of filter is chosen, the intersection of the filters (executed in the original graph) will be returned.
#' @param filterOnlyLargerComponent Optional boolean parameter that allows the preservation of the largest (with the most connections) component/group of the graph. If more than one type of filter is chosen, the intersection of the filters (executed in the original graph) will be returned.
#' @param filterDifferentNamesConnected Optional boolean parameter that allows the preservation of components (complete graphs) containing more than one species name (binomial name). If more than one type of filter is chosen, the intersection of the filters (executed in the original graph) will be returned.
#' @param filterSameNamesNotConnected Optional boolean parameter that allows the preservation of each nodes that contain the same species name (binomial name) but are not connected. If more than one type of filter is chosen, the intersection of the filters (executed in the original graph) will be returned.
#'
#' @return A list that holds relevant data of the clustering. Possible members of the list are described below.
#'
#' graph: An igraph object graph. Could be a parameter of plotc() method to plot the desired cluster(s).
#'
#' maxCliques: The largest subset of nodes in which each node is directly connected to every other node in the subset. An example would be all the possible species groups that could be delimited in the graph, which could result in groups having genomes in common.
#'
#' components: Contains the isolated nodes or groups formed of complete graphs.
#'
#' @examples
#' // Example 1
#' basicResult <- prokluster(filenames = "ANIb_percentage_identity.tab", cutoffs = 0.9)
#'
#' // Example 2
#' files <- c("ANIb_percentage_identity.tab", "ANIb_alignment_coverage.tab")
#' thresholds <- c(0.95, 0.70)
#' renamedResults1 <- prokluster(filenames = files, cutoffs = thresholds, nodesDictionary = "dictionary.tab", filterRemoveIsolated = TRUE)
#'
#' // Example 3
#' renamedResults2 <- prokluster(filenames = files, cutoffs = thresholds, nodesDictionary = "dictionary.tab", filterDifferentNamesConnected = TRUE)
#'
#' // Example 4
#' nodesNames <- read.table(file= "dictionary.tab", sep = "\t", header = F, stringsAsFactors=FALSE)
#' renamedResults3 <- prokluster(filenames = files, cutoffs = thresholds, nodesPreviousNames = nodesNames$V1, nodesTranslatedNames = nodesNames$V2, filterSameNamesNotConnected = T)
#'
#' @export
prokluster <- function(filenames,
                  cutoffs,
                  nodesDictionary = NULL,
                  nodesPreviousNames = NULL,
                  nodesTranslatedNames = NULL,
                  filterRemoveIsolated = FALSE,
                  filterRemoveLargerComponent = FALSE,
                  filterOnlyLargerComponent = FALSE,
                  filterDifferentNamesConnected = FALSE,
                  filterSameNamesNotConnected = FALSE)
{
  # Throw error if parameters have incorrect combinations
  if (!is.null(nodesDictionary) && (!is.null(nodesPreviousNames) || !is.null(nodesTranslatedNames)))
  {
    stop("If you wish to rename some nodes, please set either nodesDictionary parameter or nodesPreviousNames and nodesTranslatedNames parameters.")
  }

  if ((is.null(nodesPreviousNames) && !is.null(nodesTranslatedNames)) || (is.null(nodesTranslatedNames) && !is.null(nodesPreviousNames)))
  {
    stop("Please set both nodesPreviousNames and nodesTranslatedNames parameters.")
  }

  if (isTRUE(filterRemoveLargerComponent) && isTRUE(filterOnlyLargerComponent))
  {
    stop("Incompatible filter combination.")
  }

  finalMatrix <- matricesCombination(filenames, cutoffs)

  if (!is.null(nodesDictionary))
  {
    nodesNames <- read.table(file = nodesDictionary, sep = "\t", header = F, stringsAsFactors=FALSE)
    nodesPreviousNames <- nodesNames$V1
    nodesTranslatedNames <-  nodesNames$V2
  }

  if (!is.null(nodesPreviousNames) && !is.null(nodesTranslatedNames))
  {
    mm <- match(rownames(finalMatrix), nodesPreviousNames)
    rownames(finalMatrix)[!is.na(mm)] <- as.character(nodesTranslatedNames[na.omit(mm)])
    colnames(finalMatrix) <- rownames((finalMatrix))
  }

  library(igraph)

  # Tell Igraph it is an adjency matrix with default parameters
  g = graph_from_adjacency_matrix(as.matrix(finalMatrix))
  g <- runGraphFilters(g, filterRemoveIsolated, filterRemoveLargerComponent, filterOnlyLargerComponent, filterDifferentNamesConnected, filterSameNamesNotConnected)

  oldWarningOption <- getOption("warn")
  options(warn = -1)

  maxCliques <- max_cliques(g)

  options(warn = oldWarningOption)

  capture.output(print(maxCliques), file="maxCliques.txt")

  components <- as.data.frame(components(g)$membership)
  write.csv(components, file="components.csv")

  results <- list(graph = g, maxCliques = maxCliques, components = components)

  return(results)
}
