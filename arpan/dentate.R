graphics.off()

library(reticulate)
library(veloviz)
library(velocyto.R)
library(dplyr)

scv = import("scvelo")
adata = scv$datasets$dentategyrus()

#extract spliced, unspliced counts
spliced <- as.matrix(Matrix::t(adata$layers['spliced']))
unspliced <- as.matrix(Matrix::t(adata$layers['unspliced']))
cells <- adata$obs_names$values
genes <- adata$var_names$values
colnames(spliced) <- colnames(unspliced) <- cells
rownames(spliced) <- rownames(unspliced) <- genes

#clusters
clusters <- adata$obs$clusters
names(clusters) <- adata$obs_names$valuesnc

# density based subsampling
set.seed(0)
good.cells <- tibble::rownames_to_column(as.data.frame(clusters), "cell")
good.cells <- tibble(good.cells) %>% group_by(`clusters`) %>% sample_n(100, replace = TRUE)
good.cells <- unique(unlist(good.cells[1]))
spliced <- spliced[,good.cells]
unspliced <- unspliced[,good.cells]
clusters <- clusters[good.cells]

dim(spliced)
dim(unspliced)

#keep genes with >10 total counts
good.genes <- genes[rowSums(spliced) > 10 & rowSums(unspliced) > 10]
spliced <- spliced[good.genes,]
unspliced <- unspliced[good.genes,]

dim(spliced)
dim(unspliced)

counts <- spliced + unspliced # use combined spliced and unspliced counts
cpm <- normalizeDepth(counts) # normalize to counts per million
varnorm <- normalizeVariance(cpm) # variance stabilize, find overdispersed genes
lognorm <- log10(varnorm + 1) # log normalize

#PCA on centered and scaled expression of overdispersed genes
pcs <- reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)

#cell distance in PC space
cell.dist <- as.dist(1 - cor(t(pcs))) 


# FIXME: different params??
vel <- gene.relative.velocity.estimates(spliced,
                                       unspliced,
                                       kCells = 30,
                                       cell.dist = cell.dist,
                                       fit.quantile = 0.1)


#choose colors based on clusters for plotting later
# FIXME: how many n for rainbow ??
cell.cols <- rainbow(14)[as.numeric(clusters)]
names(cell.cols) <- names(clusters)

curr <- vel$current
proj <- vel$projected

# FIXME: different params??
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
  k = 5, # large k
  similarity.threshold = -0.5, # low t_t
  distance.weight = 0, # small omega
  distance.threshold = 1, # high d_t
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

par(mfrow = c(2,2))

#PCA
emb.pca <- pcs[,1:2]
plotEmbedding(emb.pca, colors = cell.cols, main='PCA',
              xlab = "PC1", ylab = "PC2")

#tSNE
set.seed(0)
emb.tsne <- Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) <- rownames(pcs)
plotEmbedding(emb.tsne, colors = cell.cols, main='tSNE',
              xlab = "t-SNE X", ylab = "t-SNE Y")

##UMAP
set.seed(0)
emb.umap <- uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, colors = cell.cols, main='UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

#veloviz
emb.veloviz <- veloviz$fdg_coords
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')

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

