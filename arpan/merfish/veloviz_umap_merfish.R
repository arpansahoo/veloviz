graphics.off()

library(veloviz)
library(velocyto.R)
source("../as_nn_graph.R")

# Load preprocessed data
# MERFISH data from Xia et. al., PNAS, 2019. This data is provided with the VeloViz package.
col <- MERFISH$col
pcs <- MERFISH$pcs
vel <- MERFISH$vel

curr <- vel$current
proj <- vel$projected

# Load cell cycle genes
merfish.genes <- rownames(curr)

# GO cell cycle genes (GO:0000278)
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_MITOTIC_CELL_CYCLE
cycle.genes.go <- read.csv("geneset.csv",header = FALSE)
cycle.genes.go <- cycle.genes.go$V1

merfish.cycle.go <- merfish.genes[which(merfish.genes %in% cycle.genes.go)]
curr.go <- curr[merfish.cycle.go,]
proj.go <- proj[merfish.cycle.go,]

# MERFISH genes exhibiting cell-cycle-dependent expression (Xia et al 2019, Supp Dataset 8)
# https://www.pnas.org/content/116/39/19490
cycle.genes.pnas <- read.csv("pnas.01.csv",header = TRUE)
cycle.genes.pnas <- cycle.genes.pnas$Gene

merfish.cycle.pnas <- merfish.genes[which(merfish.genes %in% cycle.genes.pnas)]
curr.pnas <- curr[merfish.cycle.pnas,]
proj.pnas <- proj[merfish.cycle.pnas,]

# build VeloViz embedding using all genes
veloviz.all <- buildVeloviz(
  curr = curr,
  proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 100,
  similarity.threshold = 0,
  distance.weight = 1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

# Plot veloviz
emb.all.vv = veloviz.all$fdg_coords
plotEmbedding(emb.all.vv, colors = col[rownames(emb.all.vv)], main = 'VeloViz')

# Normal UMAP
set.seed(0)
emb.all.umap = uwot::umap(X = pcs[,1:5], min_dist = 0.3)
rownames(emb.all.umap) = rownames(pcs)
plotEmbedding(emb.all.umap, colors = col, main = 'UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

# UMAP (initialized with velo)
set.seed(0)
nnGraph.all <- as_nn_graph(graph = veloviz.all$graph, k = 100)
emb.all.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph.all, min_dist = 0.3)
rownames(emb.all.umapVelo) <- rownames(emb.all.vv)
plotEmbedding(emb.all.umapVelo, colors = col[rownames(emb.all.umapVelo)], main = 'UMAP (initialized with VeloViz)',
              xlab = "UMAP X", ylab = "UMAP Y")

# show velocities
show.velocity.on.embedding.cor(scale(emb.all.vv), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=col[rownames(emb.all.vv)], main='VeloViz')

show.velocity.on.embedding.cor(scale(emb.all.umap), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=col, main='UMAP')

show.velocity.on.embedding.cor(scale(emb.all.umapVelo), vel,
                               n = 50,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, do.par = F,
                               cell.colors = col[rownames(emb.all.umapVelo)],
                               main='UMAP (initialized with VeloViz)')

# # Consistency scores
# deltaExp.all <- vel$deltaE
# score.all.veloviz <- consistency(emb.all.vv, deltaExp.all, nNeighbors = 10, plot.hist = FALSE)
# score.all.umapVelo <- consistency(emb.all.umapVelo, deltaExp.all, nNeighbors = 10, plot.hist = FALSE)
# score.all.umap <- consistency(emb.all.umap, deltaExp.all, nNeighbors = 10, plot.hist = FALSE)

# ks.test(score.all.veloviz, score.all.umapVelo, alternative = "two.sided") # are veloviz and umap-velo different?
# ks.test(score.all.veloviz, score.all.umapVelo, alternative = "greater") # is veloviz better than umap-velo?
# ks.test(score.all.veloviz, score.all.umap, alternative = "greater") # is veloviz better than normal UMAP?
# ks.test(score.all.umapVelo, score.all.umap, alternative = "greater") # is umap-velo better than normal UMAP?


# # Build VeloViz embedding using GO cell cycle genes
# # first, reduce dimensions.
# curr.go.norm = normalizeDepth(curr.go)
# curr.go.norm = normalizeVariance(curr.go.norm, details = TRUE)
# curr.go.norm = log10(curr.go.norm$matnorm + 1)
# curr.go.pca = RSpectra::svds(A = t(as.matrix(curr.go.norm)), k = 50,
#                              opts = list(center = TRUE, scale = FALSE,
#                                          maxitr = 2000, tol = 1e-10))
# curr.go.pca = curr.go.pca$u
# rownames(curr.go.pca) = rownames(pcs)
# 
# # Now, build embeddings.
# veloviz.go <- buildVeloviz(
#   curr = curr.go,
#   proj = proj.go,
#   normalize.depth = TRUE,
#   use.ods.genes = FALSE,
#   pca = TRUE,
#   nPCs = 3,
#   center = TRUE,
#   scale = TRUE,
#   k = 20,
#   similarity.threshold = 0,
#   distance.weight = 0.1,
#   distance.threshold = 1,
#   weighted = TRUE,
#   seed = 0,
#   verbose = FALSE
# )
# 
# # Plot veloviz
# emb.go.vv <- veloviz.go$fdg_coords
# plotEmbedding(emb.go.vv, colors = col[rownames(emb.go.vv)],
#               main = 'GO cell cycle genes - veloviz')
# 
# # Normal UMAP
# set.seed(0)
# emb.go.umap = uwot::umap(X = curr.go.pca[,1:5], min_dist = 0.3)
# rownames(emb.go.umap) = rownames(curr.go.pca)
# plotEmbedding(emb.go.umap, colors = col, main = 'GO cell cycle genes - UMAP',
#               xlab = "UMAP X", ylab = "UMAP Y")
# 
# # UMAP (initialized with velo)
# set.seed(0)
# nnGraph.go <- as_nn_graph(graph = veloviz.go$graph, k = 20)
# emb.go.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph.go, min_dist = 0.3)
# rownames(emb.go.umapVelo) <- rownames(emb.go.vv)
# plotEmbedding(emb.go.umapVelo, colors = col[rownames(emb.go.umapVelo)], main = 'GO cell cycle genes - UMAP (initialized with VeloViz)',
#               xlab = "UMAP X", ylab = "UMAP Y")
# 
# # # Consistency scores
# # deltaExp.go <- vel$deltaE[merfish.cycle.go,]
# # score.go.veloviz <- consistency(emb.go.vv, deltaExp.go, nNeighbors = 10, plot.hist = FALSE)
# # score.go.umapVelo <- consistency(emb.go.umapVelo, deltaExp.go, nNeighbors = 10, plot.hist = FALSE)
# # score.go.umap <- consistency(emb.go.umap, deltaExp.go, nNeighbors = 10, plot.hist = FALSE)
# 
# # ks.test(score.go.veloviz, score.go.umapVelo, alternative = "two.sided") # are veloviz and umap-velo different?
# # ks.test(score.go.veloviz, score.go.umapVelo, alternative = "greater") # is veloviz better than umap-velo?
# # ks.test(score.go.veloviz, score.go.umap, alternative = "greater") # is veloviz better than normal UMAP?
# # ks.test(score.go.umapVelo, score.go.umap, alternative = "greater") # is umap-velo better than normal UMAP?
# 
# 
# # build VeloViz embedding with cell-cycle dependent genes
# # first, reduce dimensions.
# curr.pnas.norm = normalizeDepth(curr.pnas)
# curr.pnas.norm = normalizeVariance(curr.pnas.norm, details = TRUE)
# curr.pnas.norm = log10(curr.pnas.norm$matnorm + 1)
# curr.pnas.pca = RSpectra::svds(A = t(as.matrix(curr.pnas.norm)), k = 50,
#                                opts = list(center = TRUE, scale = FALSE,
#                                            maxitr = 2000, tol = 1e-10))
# curr.pnas.pca = curr.pnas.pca$u
# rownames(curr.pnas.pca) = rownames(pcs)
# 
# veloviz.pnas <- buildVeloviz(
#   curr = curr.pnas,
#   proj = proj.pnas,
#   normalize.depth = TRUE,
#   use.ods.genes = FALSE,
#   pca = TRUE,
#   nPCs = 3,
#   center = TRUE,
#   scale = TRUE,
#   k = 50,
#   similarity.threshold = 0.5,
#   distance.weight = 0.01,
#   distance.threshold = 1,
#   weighted = TRUE,
#   seed = 0,
#   verbose = FALSE
# )
# 
# # Plot veloviz 
# emb.pnas.vv <- veloviz.pnas$fdg_coords
# plotEmbedding(emb.pnas.vv, colors = col[rownames(emb.pnas.vv)],
#               main = 'Xia et al cell cycle genes - veloviz')
# 
# # Normal UMAP 
# set.seed(0)
# emb.pnas.umap = uwot::umap(X = curr.pnas.pca[,1:5], min_dist = 0.3)
# rownames(emb.pnas.umap) = rownames(curr.pnas.pca)
# plotEmbedding(emb.pnas.umap, colors = col, main = 'Xia et al cell cycle genes - UMAP',
#               xlab = "UMAP X", ylab = "UMAP Y")
# 
# # UMAP (initialized with velo)
# set.seed(0)
# nnGraph.pnas <- as_nn_graph(graph = veloviz.pnas$graph, k = 50)
# emb.pnas.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph.pnas, min_dist = 0.3)
# rownames(emb.pnas.umapVelo) <- rownames(emb.pnas.vv)
# plotEmbedding(emb.pnas.umapVelo, colors = col[rownames(emb.pnas.umapVelo)], main = 'Xia et al cell cycle genes - UMAP (initialized with VeloViz)',
#               xlab = "UMAP X", ylab = "UMAP Y")
# 
# # # Consistency scores
# # library(RANN)
# 
# # deltaExp.pnas <- vel$deltaE[merfish.cycle.pnas,]
# # score.pnas.veloviz <- consistency(emb.pnas.vv, deltaExp.pnas, nNeighbors = 10, plot.hist = FALSE)
# # score.pnas.umapVelo <- consistency(emb.pnas.umapVelo, deltaExp.pnas, nNeighbors = 10, plot.hist = FALSE)
# # score.pnas.umap <- consistency(emb.pnas.umap, deltaExp.pnas, nNeighbors = 10, plot.hist = FALSE)
# 
# # ks.test(score.pnas.veloviz, score.pnas.umapVelo, alternative = "two.sided") # are veloviz and umap-velo different?
# # ks.test(score.pnas.veloviz, score.pnas.umapVelo, alternative = "greater") # is veloviz better than umap-velo?
# # ks.test(score.pnas.veloviz, score.pnas.umap, alternative = "greater") # is veloviz better than normal UMAP?
# # ks.test(score.pnas.umapVelo, score.pnas.umap, alternative = "greater") # is umap-velo better than normal UMAP?
