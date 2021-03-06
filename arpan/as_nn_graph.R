library(igraph)

as_nn_graph = function(graph, k) {
  edges <- get.edgelist(graph) # array of [source, target] 
  weights <- E(graph)$weight # edge weights
  numEdges <- gsize(graph)
  numVertices <- gorder(graph)
  
  # array of [X, Neighbors(X)]
  idx <- matrix(nrow=numVertices, ncol=k+1) 
  
  # each row is distances from a vertex X to its neighbors
  dist <- matrix(0, nrow=numVertices, ncol=k+1)
  
  # replace NA values by saying a vertex X is a neighbor of itself 
  for (i in 1:numVertices) {
    idx[i,] = i
  }
  
  for (i in 1:numEdges) {
    edge <- edges[i,] # edge = {source, target}
    source <- edge[1]
    target <- edge[2]
    
    # index of first "NA" value in row, after index 1
    index <- min(which(idx[source,] == source)[-1])
    
    idx[source, index] = target
    dist[source, index] = weights[i]
  }
  
  # construct & return a list consisting of two elements: "idx", "dist"
  nnGraph = list()
  nnGraph$idx = idx
  nnGraph$dist = dist
  
  return(nnGraph)
}