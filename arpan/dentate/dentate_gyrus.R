library(veloviz)
library(reticulate)
library(velocyto.R)
library(viridis)
source("../as_nn_graph.R")

scv <- import("scvelo")
adata <- scv$datasets$dentategyrus()

#extract spliced, unspliced counts
spliced <- as.matrix(Matrix::t(adata$layers['spliced']))
unspliced <- as.matrix(Matrix::t(adata$layers['unspliced']))
cells <- adata$obs_names$values
genes <- adata$var_names$values
colnames(spliced) <- colnames(unspliced) <- cells
rownames(spliced) <- rownames(unspliced) <- genes

#clusters
clusters <- adata$obs$clusters
names(clusters) <- adata$obs_names$values
col = rev(plasma(length(levels(clusters))))
col.clust = col[clusters] 
names(col.clust) <- cells

pdf("dentate_legend.pdf")
par(mfrow=c(1,1))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='')
# uniqueCols <- unique(col.clust[rownames(emb.veloviz1)])
# order <- order(uniqueCols)
# order <- c(14,2,7,12,4,13,5,3,1,9,10,8,6,11)
legend("topleft", legend = unique(clusters)[order], col = uniqueCols[order], pch=16, cex=0.7, ncol=1)
dev.off()

#keep genes with >10 total counts
good.genes = genes[rowSums(spliced) > 10 & rowSums(unspliced) > 10]
spliced = spliced[good.genes,]
unspliced = unspliced[good.genes,]

#normalize
counts = spliced + unspliced # use combined spliced and unspliced counts
cpm = normalizeDepth(counts) # normalize to counts per million
varnorm = normalizeVariance(cpm) # variance stabilize, find overdispersed genes
lognorm = log10(varnorm + 1) # log normalize


#pca
pcs = reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)
cell.dist = as.dist(1-cor(t(pcs))) 

vel <- readRDS(file = "dentate_gyrus_vel01.rds")
curr <- vel$current
proj <- vel$projected


##build veloviz 

##final 
veloviz1 <- buildVeloviz(
  curr = curr, 
  proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 10,
  center = TRUE,
  scale = TRUE,
  k = 100, #100
  similarity.threshold = 0,
  distance.weight = 1,
  distance.threshold = 0.9, #0.7
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.veloviz1 = veloviz1$fdg_coords
par(mfrow=c(1,1))
plotEmbedding(emb.veloviz1, colors=col.clust[rownames(emb.veloviz1)], 
              main="VeloViz", xlab="VeloViz X", ylab="VeloViz Y")
legend(x=9.5, y=-8, legend = unique(clusters), col = unique(col.clust), pch=16, cex=0.8, ncol=2)

##UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, colors=col.clust,
              main='UMAP', xlab = "UMAP X", ylab = "UMAP Y")
legend(x=-18, y=19.5, legend = unique(clusters), col = unique(col.clust), pch=16, cex=0.7, ncol=1)



# Convert veloviz$graph (igraph type) to an idx & dist representation
nnGraph <- as_nn_graph(graph = veloviz1$graph, k = 100)

# input nnGraph to UMAP and plot
set.seed(0)
emb.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
rownames(emb.umapVelo) <- rownames(emb.veloviz1)
plotEmbedding(emb.umapVelo, colors=col.clust[rownames(emb.umapVelo)], 
              main = 'UMAP with VeloViz', xlab = "UMAP X", ylab = "UMAP Y")
legend(x=-7.5, y=4.5, legend = unique(clusters), col = unique(col.clust[rownames(emb.umapVelo)]), pch=16, cex=0.7, ncol=1)


# VELOCITIES
par(mfrow=c(2,2), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))

show.velocity.on.embedding.cor(emb.veloviz1, vel,
                               n = 50,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F, 
                               frame.plot = TRUE, xaxt='n',yaxt='n',xlab="VeloViz X",ylab="VeloViz Y",
                               cell.colors=scales::alpha(col.clust[rownames(emb.veloviz1)],0.4), main='VeloViz')

show.velocity.on.embedding.cor(emb.umap, vel,
                               n = 50,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               frame.plot = TRUE, xaxt='n',yaxt='n',xlab="UMAP X",ylab="UMAP Y",
                               cell.colors=scales::alpha(col.clust,0.4), main='UMAP')
legend(x=-18, y=19.5, legend = unique(clusters), col = unique(col.clust), pch=16, cex=0.7, ncol=1)

show.velocity.on.embedding.cor(emb.umapVelo, vel,
                               n = 50,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               frame.plot = TRUE, xaxt='n',yaxt='n',xlab="UMAP X",ylab="UMAP Y",
                               cell.colors=scales::alpha(col.clust[rownames(emb.umapVelo)], 0.4), main='UMAP with VeloViz')
legend(x=-7.5, y=4.5, legend = unique(clusters), col = unique(col.clust[rownames(emb.umapVelo)]), pch=16, cex=0.7, ncol=1)

