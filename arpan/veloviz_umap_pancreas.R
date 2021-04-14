graphics.off()

library(reticulate)
library(veloviz)
library(velocyto.R)
library(viridis)
library(igraph)

# clusters <- pancreas$clusters # cell type annotations
# pcs <- pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
# vel <- pancreas$vel # velocity

# getting pancreas data from scVelo
scv <- import("scvelo")
adata <- scv$datasets$pancreas()

# extract spliced, unspliced counts
spliced <- as.matrix(Matrix::t(adata$layers['spliced']))
unspliced <- as.matrix(Matrix::t(adata$layers['unspliced']))
cells <- adata$obs_names$values
genes <- adata$var_names$values
colnames(spliced) <- colnames(unspliced) <- cells
rownames(spliced) <- rownames(unspliced) <- genes

# clusters
clusters <- adata$obs$clusters
names(clusters) <- adata$obs_names$values

# subsample to make things faster
set.seed(0)
good.cells <- sample(cells, length(cells) / 2)
spliced <- spliced[,good.cells]
unspliced <- unspliced[,good.cells]
clusters <- clusters[good.cells]

# keep genes with >10 total counts
good.genes <- genes[rowSums(spliced) > 10 & rowSums(unspliced) > 10]
spliced <- spliced[good.genes,]
unspliced = unspliced[good.genes,]

counts <- spliced + unspliced # use combined spliced and unspliced counts
cpm <- normalizeDepth(counts) # normalize to counts per million
varnorm <- normalizeVariance(cpm) # variance stabilize, find overdispersed genes
lognorm <- log10(varnorm + 1) # log normalize

# PCA on centered and scaled expression of overdispersed genes
pcs <- reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)

# # cell distance in PC space
# cell.dist <- as.dist(1-cor(t(pcs)))
#
# # compute velocity
# vel <- gene.relative.velocity.estimates(spliced,
#                                         unspliced,
#                                         kCells = 30,
#                                         cell.dist = cell.dist,
#                                         fit.quantile = 0.1)
# saveRDS(vel, file = "pancreas_vel.rds")
vel <- readRDS(file = "pancreas_vel.rds")

# choose colors based on clusters for plotting later
col = rev(plasma(length(levels(clusters))))
cell.cols = col[clusters] 
names(cell.cols) = names(clusters)

# pdf("pancreas_legend.pdf")
# par(mfrow=c(1,1))
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='')
# uniqueCols <- unique(cell.cols)
# orderCols <- rev(order(uniqueCols))
# legend("topleft", 
#        legend = unique(clusters)[c(2, 3, 5, 1, 4, 6, 8, 7)],
#        col = uniqueCols[c(2, 3, 5, 1, 4, 6, 8, 7)],
#        pch=16, cex=0.7, ncol=1)
# dev.off()

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
  k = 20,
  similarity.threshold = 0.25,
  distance.weight = 1,
  distance.threshold = 0.5,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

par(mfrow=c(1,1), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))

# Plot veloviz
emb.veloviz = veloviz$fdg_coords
# plotEmbedding(emb.veloviz, 
#               colors=cell.cols[rownames(emb.veloviz)],
#               frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
#               main='VeloViz', xlab="VeloViz X", ylab = "VeloViz Y")
# 
# # UMAP (normal)
# set.seed(0)
# emb.umap = uwot::umap(pcs, min_dist = 0.5)
# rownames(emb.umap) <- rownames(pcs)
# plotEmbedding(emb.umap, colors = cell.cols, main='UMAP',
#               frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
#               xlab = "UMAP X", ylab = "UMAP Y")

# Convert veloviz$graph (igraph type) to an idx & dist representation
nnGraph <- veloviz::asNNGraph(veloviz)

# input nnGraph to UMAP and plot
set.seed(0)
emb.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
rownames(emb.umapVelo) <- rownames(emb.veloviz)
plotEmbedding(emb.umapVelo,
              colors = cell.cols[rownames(emb.umapVelo)],
              frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
              main = 'UMAP with VeloViz', xlab = "UMAP X", ylab = "UMAP Y")

pdf("pancreas_new.pdf")

# show velocities
# par(mfrow=c(2,2), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))
# show.velocity.on.embedding.cor(scale(emb.veloviz), vel,
#                                n = 50,
#                                scale='sqrt',
#                                cex=1, arrow.scale=1, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
#                                frame.plot = TRUE, xaxt = 'n', yaxt = 'n', xlab="VeloViz X", ylab='VeloViz Y',
#                                cell.colors=scales::alpha(cell.cols[rownames(emb.veloviz)], 0.4), main='VeloViz')
# show.velocity.on.embedding.cor(scale(emb.umap), vel,
#                                n = 50,
#                                scale='sqrt',
#                                cex=1, arrow.scale=1, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
#                                frame.plot = TRUE, xaxt = 'n', yaxt = 'n', xlab="UMAP X", ylab='UMAP Y',
#                                cell.colors=scales::alpha(cell.cols, 0.4), main='UMAP')
show.velocity.on.embedding.cor(scale(emb.umapVelo), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1.5, do.par = F,
                               cell.colors = scales::alpha(cell.cols[rownames(emb.umapVelo)], 0.4),
                               frame.plot = TRUE, xaxt = 'n', yaxt = 'n', xlab="UMAP X", ylab='UMAP Y',
                               main='UMAP with VeloViz')

dev.off()
