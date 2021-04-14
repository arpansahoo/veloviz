library(igraph)

asNNGraph = function(vig) {
  graph <- vig$graph
  edges <- igraph::get.edgelist(graph) 
  numEdges <- igraph::gsize(graph)
  numVertices <- igraph::gorder(graph)
  
  k <- max(igraph::degree(graph, mode="out")) # k is the max out-degree
  allDists <- vig$projected_neighbors$all_dists # distances
  
  # each row i contains the indices of vertex i's NNs
  # vertex i is a NN of itself so first value in row is i
  idx <- matrix(nrow=numVertices, ncol=k+1) 
  
  # each row i contains distances from vertex i to its NNs
  # vertex i is a NN of itself so first value in row is 0.0
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
  
  # return a list consisting of two elements: "idx" and "dist"
  out <- list()
  out[['idx']] <- idx
  out[['dist']] <- dist  
  return(out)
}