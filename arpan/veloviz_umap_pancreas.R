graphics.off()

library(reticulate)
library(veloviz)
library(velocyto.R)
library(tictoc)
source("as_nn_graph.R")

# clusters <- pancreas$clusters # cell type annotations
# pcs <- pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
# vel <- pancreas$vel # velocity

#getting pancreas data from scVelo
scv <- import("scvelo")
adata <- scv$datasets$pancreas()

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

#subsample to make things faster
set.seed(0)
good.cells <- sample(cells, length(cells) / 2)
spliced <- spliced[,good.cells]
unspliced <- unspliced[,good.cells]
clusters <- clusters[good.cells]

#keep genes with >10 total counts
good.genes <- genes[rowSums(spliced) > 10 & rowSums(unspliced) > 10]
spliced <- spliced[good.genes,]
unspliced = unspliced[good.genes,]

counts <- spliced + unspliced # use combined spliced and unspliced counts
cpm <- normalizeDepth(counts) # normalize to counts per million
varnorm <- normalizeVariance(cpm) # variance stabilize, find overdispersed genes
lognorm <- log10(varnorm + 1) # log normalize

#PCA on centered and scaled expression of overdispersed genes
pcs <- reduceDimensions(lognorm, center = TRUE, scale = TRUE, nPCs = 50)

#cell distance in PC space
cell.dist <- as.dist(1-cor(t(pcs))) # cell distance in PC space

# compute velocity
vel <- gene.relative.velocity.estimates(spliced,
                                        unspliced,
                                        kCells = 30,
                                        cell.dist = cell.dist,
                                        fit.quantile = 0.1)


# choose colors based on clusters for plotting later
cell.cols <- rainbow(8)[as.numeric(clusters)]
names(cell.cols) <- names(clusters)

# build VeloViz embedding
curr <- vel$current
proj <- vel$projected
# veloviz <- buildVeloviz(
#   curr = curr, 
#   proj = proj,
#   normalize.depth = TRUE,
#   use.ods.genes = TRUE,
#   alpha = 0.05,
#   pca = TRUE,
#   nPCs = 20,
#   center = TRUE,
#   scale = TRUE,
#   k = 5,
#   similarity.threshold = 0.25,
#   distance.weight = 1,
#   distance.threshold = 0.5,
#   weighted = TRUE,
#   seed = 0,
#   verbose = FALSE
# )
# emb.veloviz = veloviz$fdg_coords

#normalize depth
curr.norm = normalizeDepth(curr)
proj.norm = normalizeDepth(proj)

#variance stabilize current
curr.varnorm.info = normalizeVariance(curr.norm, details = TRUE)
curr.varnorm = curr.varnorm.info$matnorm

#use same model for projected
scale.factor = curr.varnorm.info$df$scale_factor #gene scale factors
names(scale.factor) = rownames(curr.varnorm.info$df)

m = proj.norm
rmean = Matrix::rowMeans(m) #row mean
sumx = Matrix::rowSums(m)
sumxx = Matrix::rowSums(m^2)
rsd = sqrt((sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2) / (ncol(m)-1)) #row sd

proj.varnorm = proj.norm / rsd * scale.factor[names(rsd)]
proj.varnorm = proj.norm[rownames(curr.varnorm),]

#log normalize
curr.pca = log10(curr.varnorm + 1)
proj.pca = log10(proj.varnorm + 1)

#mean center
c.rmean = Matrix::rowMeans(curr.pca)
curr.pca = curr.pca - c.rmean
p.rmean = Matrix::rowMeans(proj.pca)
proj.pca = proj.pca - p.rmean

#scale variance
c.sumx = Matrix::rowSums(curr.pca)
c.sumxx = Matrix::rowSums(curr.pca^2)
c.rsd = sqrt((c.sumxx - 2*c.sumx*c.rmean + ncol(curr.pca)*c.rmean^2)/(ncol(curr.pca)-1))
curr.pca = curr.pca/c.rsd

p.sumx = Matrix::rowSums(proj.pca)
p.sumxx = Matrix::rowSums(proj.pca^2)
p.rsd = sqrt((p.sumxx - 2*p.sumx*p.rmean + ncol(proj.pca)*p.rmean^2)/(ncol(proj.pca)-1))
proj.pca = proj.pca/p.rsd

#PCA
pca = RSpectra::svds(A = Matrix::t(curr.pca), k=20,
                     opts = list(
                       center = FALSE, ## already done
                       scale = FALSE, ## already done
                       maxitr = 2000,
                       tol = 1e-10))

#scores of current and projected
curr.scores = Matrix::t(curr.pca) %*% pca$v[,1:10]
proj.scores = Matrix::t(proj.pca) %*% pca$v[,1:10]

#VeloViz graph parameters
k = 5
similarity.threshold = 0.25
distance.weight = 1
distance.threshold = 0.5
weighted = TRUE

#build graph
tic('Veloviz')
set.seed(0)
veloviz = graphViz(t(curr.scores), t(proj.scores), k,
                   cell.colors=NA,
                   similarity_threshold=similarity.threshold,
                   distance_weight = distance.weight,
                   distance_threshold = distance.threshold,
                   weighted = weighted,
                   plot = FALSE,
                   return_graph = TRUE)

emb.veloviz = veloviz$fdg_coords
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')
toc()

# Convert veloviz$graph (igraph type) to an idx & dist representation
nnGraph <- as_nn_graph(graph = veloviz$graph, k = 5)

# input nnGraph to UMAP and plot
tic('UMAP')
par(mfrow = c(1,1))
set.seed(0)
emb.umap <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
rownames(emb.umap) <- rownames(emb.veloviz)
plotEmbedding(emb.umap, colors = cell.cols[rownames(emb.umap)], main='UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")
toc()

# # alternatively, we can send in pcs as input,
# # and send nnGraph as target data for supervised dimension reduction
# emb.umap <- uwot::umap(X = pcs, y = nnGraph, min_dist = 0.5)
# 
# # show velocities
# par(mfrow = c(1,1))
# show.velocity.on.embedding.cor(scale(emb.umap), vel,
#                                n = 50,
#                                scale='sqrt',
#                                cex=1, arrow.scale=1, show.grid.flow=TRUE,
#                                min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
#                                cell.colors=cell.cols[rownames(emb.umap)], main='UMAP')
