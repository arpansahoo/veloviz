library(veloviz)
library(reticulate)
library(velocyto.R)

scv <- import("scvelo")
plt  = import('matplotlib.pyplot')
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
col.clust <- rainbow(length(table(clusters)))[as.numeric(clusters)]
names(col.clust) <- cells


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


par(mfrow=c(2,2))
emb.pca <- pcs[,1:2]
plot(emb.pca, pch=16, col=col.clust)

#tSNE
set.seed(0)
emb.tsne = Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) = rownames(pcs)
plot(emb.tsne, main='tSNE',pch = 16, col=col.clust,
     xlab = "t-SNE X", ylab = "t-SNE Y")

##UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plot(emb.umap, main='UMAP',pch = 16, col=col.clust,
     xlab = "UMAP X", ylab = "UMAP Y")

##diffusion map 
set.seed(0)
diffmap = destiny::DiffusionMap(pcs)
emb.diffmap = destiny::eigenvectors(diffmap)[,1:2]
row.names(emb.diffmap) = cells
plot(emb.diffmap, pch = 16, col=col.clust, main = "Diffusion map")

##velocity 
# vel = gene.relative.velocity.estimates(spliced,
#                                        unspliced,
#                                        kCells = 30,
#                                        cell.dist = cell.dist,
#                                        fit.quantile = 0.1)
# saveRDS(vel, file = "dentate_gyrus_vel01.rds")
vel <- readRDS(file = "dentate_gyrus_vel01.rds")
curr <- vel$current
proj <- vel$projected


##build veloviz 

##final 
veloviz1 <- buildVeloviz(
  curr = curr, proj = proj,
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
plot(emb.veloviz1, pch=16, col=col.clust[rownames(emb.veloviz1)])
legend(x=9.5, y=-7.5, legend = unique(clusters), col = unique(col.clust), pch=16, cex=0.8, ncol=3)

show.velocity.on.embedding.cor(emb.veloviz1, vel, 
                               n = 50,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')
legend(x=10, y=-8, legend = unique(clusters), col = unique(col.clust), pch=16, cex=0.7, ncol=2)


