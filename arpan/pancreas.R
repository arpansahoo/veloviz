graphics.off()

library(veloviz)
library(velocyto.R)

clusters = pancreas$clusters # cell type annotations
pcs = pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
vel = pancreas$vel # velocity

# choose colors based on clusters for plotting later
cell.cols = rainbow(8)[as.numeric(clusters)]
names(cell.cols) = names(clusters)

# buildVeloviz
curr = vel$current
proj = vel$projected

veloviz = buildVeloviz(
  curr = curr, proj = proj,
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

# Compare embeddings

par(mfrow = c(2,2))

#PCA
emb.pca = pcs[,1:2]
plotEmbedding(emb.pca, colors = cell.cols, main='PCA',
              xlab = "PC1", ylab = "PC2")

#tSNE
set.seed(0)
emb.tsne = Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) = rownames(pcs)
plotEmbedding(emb.tsne, colors = cell.cols, main='tSNE',
              xlab = "t-SNE X", ylab = "t-SNE Y")

#UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, colors = cell.cols, main='UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

#veloviz
emb.veloviz = veloviz$fdg_coords
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='veloviz')


# Project velocity using velocyto
par(mfrow = c(2,2))

show.velocity.on.embedding.cor(scale(emb.pca), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par = F,
                               cell.colors=cell.cols, main='PCA')
show.velocity.on.embedding.cor(scale(emb.tsne), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='tSNE')
show.velocity.on.embedding.cor(scale(emb.umap), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='UMAP')
show.velocity.on.embedding.cor(scale(emb.veloviz), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=cell.cols, main='VeloViz')

