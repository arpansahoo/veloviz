library(igraph)

as_nn_graph = function(veloviz) {
  graph <- veloviz$graph
  edges <- get.edgelist(graph) # array of [source, target] 
  numEdges <- gsize(graph)
  numVertices <- gorder(graph)
  
  allDists <- veloviz$projected_neighbors$all_dists # distances
  k <- max(degree(veloviz$graph, mode="out"))
  
  # array of [X, Neighbors(X)]
  idx <- matrix(nrow=numVertices, ncol=k+1) 
  
  # each row is distances from a vertex X to its neighbors
  dist <- matrix(0, nrow=numVertices, ncol=k+1)
  
  # replace NA values by saying a vertex X is a neighbor of itself 
  for (i in 1:numVertices) {
    idx[i,] <- i
  }
  
  for (i in 1:numEdges) {
    edge <- edges[i,] # edge = {source, target}
    source <- edge[1]
    target <- edge[2]
    
    # index of first "NA" value in row, after index 1
    index <- min(which(idx[source,] == source)[-1])
    
    idx[source, index] <- target
    dist[source, index] <- allDists[source, target]
  }
  
  # construct & return a list consisting of two elements: "idx", "dist"
  nnGraph <- list()
  nnGraph$idx <- idx
  nnGraph$dist <- dist
  
  return(nnGraph)
}