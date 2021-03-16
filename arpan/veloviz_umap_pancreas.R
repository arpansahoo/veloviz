graphics.off()

library(reticulate)
library(veloviz)
library(RANN)
library(velocyto.R)
library(tictoc)
source("as_nn_graph.R")

clusters <- pancreas$clusters # cell type annotations
pcs <- pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
vel <- pancreas$vel # velocity

# # getting pancreas data from scVelo
# scv <- import("scvelo")
# adata <- scv$datasets$pancreas()
# 
# # extract spliced, unspliced counts
# spliced <- as.matrix(Matrix::t(adata$layers['spliced']))
# unspliced <- as.matrix(Matrix::t(adata$layers['unspliced']))
# cells <- adata$obs_names$values
# genes <- adata$var_names$values
# colnames(spliced) <- colnames(unspliced) <- cells
# rownames(spliced) <- rownames(unspliced) <- genes
# 
# # clusters
# clusters <- adata$obs$clusters
# names(clusters) <- adata$obs_names$values
# 
# # subsample to make things faster
# set.seed(0)
# good.cells <- sample(cells, length(cells) / 2)
# spliced <- spliced[,good.cells]
# unspliced <- unspliced[,good.cells]
# clusters <- clusters[good.cells]
# 
# # keep genes with >10 total counts
# good.genes <- genes[rowSums(spliced) > 10 & rowSums(unspliced) > 10]
# spliced <- spliced[good.genes,]
# unspliced = unspliced[good.genes,]
# 
# counts <- spliced + unspliced # use combined spliced and unspliced counts
# cpm <- normalizeDepth(counts) # normalize to counts per million
# varnorm <- normalizeVariance(cpm) # variance stabilize, find overdispersed genes
# lognorm <- log10(varnorm + 1) # log normalize
# 
# # PCA on centered and scaled expression of overdispersed genes
# pcs <- reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)
# 
# # cell distance in PC space
# cell.dist <- as.dist(1-cor(t(pcs)))
# 
# # compute velocity
# vel <- gene.relative.velocity.estimates(spliced,
#                                         unspliced,
#                                         kCells = 30,
#                                         cell.dist = cell.dist,
#                                         fit.quantile = 0.1)

# choose colors based on clusters for plotting later
cell.cols <- rainbow(8)[as.numeric(clusters)]
names(cell.cols) <- names(clusters)

# build VeloViz embedding
curr <- vel$current
proj <- vel$projected
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
  k = 5,
  similarity.threshold = 0.25,
  distance.weight = 1,
  distance.threshold = 0.5,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

# Plot veloviz
emb.veloviz = veloviz$fdg_coords
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')

# UMAP (normal)
set.seed(0)
emb.normalUMAP = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.normalUMAP) <- rownames(pcs)
plotEmbedding(emb.normalUMAP, colors = cell.cols, main='UMAP (normal)',
              xlab = "UMAP X", ylab = "UMAP Y")

# Convert veloviz$graph (igraph type) to an idx & dist representation
nnGraph <- as_nn_graph(graph = veloviz$graph, k = 5)

# input nnGraph to UMAP and plot
set.seed(0)
emb.umap <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
rownames(emb.umap) <- rownames(emb.veloviz)
plotEmbedding(emb.umap,
              colors = cell.cols[rownames(emb.umap)],
              main = 'UMAP (initialized with veloviz)', xlab = "UMAP X", ylab = "UMAP Y")

# show velocities
show.velocity.on.embedding.cor(scale(emb.veloviz), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='VeloViz')

show.velocity.on.embedding.cor(scale(emb.normalUMAP), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='UMAP (normal)')

show.velocity.on.embedding.cor(scale(emb.umap), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par = F,
                               cell.colors = cell.cols[rownames(emb.umap)],
                               main='UMAP (initialized with veloviz)')

# Consistency scores
score.veloviz <- consistency(emb.veloviz, vel$deltaE, nNeighbors = 10, plot.hist = TRUE)
score.umap <- consistency(emb.umap, vel$deltaE, nNeighbors = 10, plot.hist = TRUE)
score.normalUMAP <- consistency(emb.normalUMAP, vel$deltaE, nNeighbors = 10, plot.hist = TRUE)

ks.test(score.veloviz, score.umap, alternative = "two.sided") # are veloviz and umap-velo different? 
ks.test(score.veloviz, score.normalUMAP, alternative = "greater") # is veloviz better than normal UMAP?
ks.test(score.umap, score.normalUMAP, alternative = "greater") # is umap-velo better than normal UMAP?
