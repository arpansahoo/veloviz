graphics.off()

library(veloviz)
library(velocyto.R)
library(igraph)
source("as_nn_graph.R")

clusters <- pancreas$clusters # cell type annotations
pcs <- pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
vel <- pancreas$vel # velocity

# choose colors based on clusters for plotting later
cell.cols <- rainbow(8)[as.numeric(clusters)]
names(cell.cols) <- names(clusters)

# buildVeloviz
curr <- vel$current
proj <- vel$projected
k <- 5

veloviz <- buildVeloviz(
  curr = curr, 
  proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = TRUE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 20,
  center = TRUE,
  scale = TRUE,
  k = k,
  similarity.threshold = 0.25,
  distance.weight = 1,
  distance.threshold = 0.5,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

# Convert veloviz$graph (igraph type) to an idx & dist representation
nnGraph <- as_nn_graph(graph = veloviz$graph, k = k)

# input nnGraph to UMAP and plot
par(mfrow = c(1,1))
set.seed(0)
emb.umap <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, colors = cell.cols, main='UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

# show velocities 
par(mfrow = c(1,1))
show.velocity.on.embedding.cor(scale(emb.umap), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='UMAP')