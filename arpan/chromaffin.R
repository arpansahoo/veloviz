graphics.off()

library(veloviz)
library(velocyto.R)

# Read resulting RDS file from reading in count matrices from the loom file
ldat <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/ldat.rds"))

# Reduce the cell names to the short well labels
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("_unique.bam", "", gsub(".*:", "", colnames(x)))
  x
})

# Read in cell cluster assignment used in Furlan et al. (Science’17)
cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))

# Exonic read (spliced) expression matrix
spliced <- ldat$spliced
# Intronic read (unspliced) expression matrix
unspliced <- ldat$unspliced

# Filter expression matrices based on some minimum max-cluster averages
spliced <- filter.genes.by.cluster.expression(spliced,
                                              cell.colors,
                                              min.max.cluster.average = 5)
unspliced <- filter.genes.by.cluster.expression(unspliced,
                                                cell.colors,
                                                min.max.cluster.average = 1)

goodGenes <- intersect(rownames(spliced), rownames(unspliced))
spliced <- spliced[goodGenes,]
unspliced <- unspliced[goodGenes,]

# Normalize matrix
counts <- spliced + unspliced # use combined spliced and unspliced counts
cpm <- normalizeDepth(counts) # normalize to counts per million
varnorm <- normalizeVariance(cpm) # variance stabilize, find overdispersed genes
lognorm <- log10(varnorm + 1) # log normalize

# PCA on centered and scaled expression of overdispersed genes
pcs <- reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)

# Cell distance in PC space
cell.dist <- as.dist(1 - cor(t(pcs))) 

# Compute velocity
# FIXME: different params??
vel <- gene.relative.velocity.estimates(spliced,
                                        unspliced,
                                        kCells = 30,
                                        cell.dist = cell.dist,
                                        fit.quantile = 0.1)

# Run Veloviz
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
  k = 10, #large k 
  similarity.threshold = 0.25, #low t_t
  distance.weight = 1, # small omega
  distance.threshold = 0.25, # high d_t
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

par(mfrow = c(2,2))

#PCA
emb.pca <- pcs[,1:2]
plotEmbedding(emb.pca, colors = cell.colors, main= 'PCA',
              xlab = "PC1", ylab = "PC2")

#tSNE
set.seed(0)
emb.tsne <- Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) <- rownames(pcs)
plotEmbedding(emb.tsne, colors = cell.colors, main = 'tSNE',
              xlab = "t-SNE X", ylab = "t-SNE Y")

#UMAP
set.seed(0)
emb.umap <- uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, colors = cell.colors, main = 'UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

# Veloviz
emb.veloviz <- veloviz$fdg_coords
plotEmbedding(emb.veloviz, 
              colors = cell.colors[rownames(emb.veloviz)], 
              main='veloviz')

# Project velocity using velocyto
par(mfrow = c(2,2))

show.velocity.on.embedding.cor(scale(emb.pca), vel,
                               n=50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par=F,
                               cell.colors=cell.colors, main='PCA')
show.velocity.on.embedding.cor(scale(emb.tsne), vel,
                               n=50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par=F,
                               cell.colors=cell.colors, main='tSNE')
show.velocity.on.embedding.cor(scale(emb.umap), vel,
                               n=50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par=F,
                               cell.colors=cell.colors, main='UMAP')

show.velocity.on.embedding.cor(scale(emb.veloviz), vel,
                               n=50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par=F,
                               cell.colors=cell.colors, main='VeloViz')
