graphics.off()

library(veloviz)
source("as_nn_graph.R")

# MERFISH data from Xia et. al., PNAS, 2019. This data is provided with the VeloViz package.
col = MERFISH$col
pcs = MERFISH$pcs
vel = MERFISH$vel

curr = vel$current
proj = vel$projected

# Load cell cycle genes
merfish.genes = rownames(curr)

# GO cell cycle genes (GO:0000278)
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_MITOTIC_CELL_CYCLE
cycle.genes.go = read.csv("geneset.csv",header = FALSE)
cycle.genes.go = cycle.genes.go$V1

merfish.cycle.go =  merfish.genes[which(merfish.genes %in% cycle.genes.go)]
curr.go = curr[merfish.cycle.go,]
proj.go = proj[merfish.cycle.go,]

# MERFISH genes exhibiting cell-cycle-dependent expression (Xia et al 2019, Supp Dataset 8)
# https://www.pnas.org/content/116/39/19490
cycle.genes.pnas = read.csv("pnas.01.csv",header = TRUE)
cycle.genes.pnas = cycle.genes.pnas$Gene

merfish.cycle.pnas = merfish.genes[which(merfish.genes %in% cycle.genes.pnas)]
curr.pnas = curr[merfish.cycle.pnas,]
proj.pnas = proj[merfish.cycle.pnas,]

par(mfrow = c(1,1))

# build VeloViz embedding using all genes
veloviz.all = buildVeloviz(
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

# UMAP
set.seed(0)
nnGraph.all <- as_nn_graph(graph = veloviz.all$graph, k = 100)

emb.all.umap <- uwot::umap(X = NULL, nn_method = nnGraph.all, min_dist = 0.3)
rownames(emb.all.umap) = rownames(pcs)
plotEmbedding(emb.all.umap, colors = col, main = 'all genes - UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

# build VeloViz embedding using GO cell cycle genes
curr.go.norm = normalizeDepth(curr.go)
curr.go.norm = normalizeVariance(curr.go.norm, details = TRUE)

curr.go.norm = log10(curr.go.norm$matnorm + 1)
curr.go.pca = RSpectra::svds(A = t(as.matrix(curr.go.norm)), k = 50,
                             opts = list(center = TRUE, scale = FALSE,
                                         maxitr = 2000, tol = 1e-10))
curr.go.pca = curr.go.pca$u
rownames(curr.go.pca) = rownames(pcs)

veloviz.go = buildVeloviz(
  curr = curr.go,
  proj = proj.go,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 20,
  similarity.threshold = 0,
  distance.weight = 0.1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

# UMAP
par(mfrow = c(1,1))
set.seed(0)
nnGraph.go <- as_nn_graph(graph = veloviz.go$graph, k = 20)

emb.go.umap <- uwot::umap(X = NULL, nn_method = nnGraph.go, min_dist = 0.3)
rownames(emb.go.umap) = rownames(curr.go.pca)
plotEmbedding(emb.go.umap, colors = col, main = 'GO cell cycle genes - UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

# build VeloViz embedding with cell-cycle dependent genes
curr.pnas.norm = normalizeDepth(curr.pnas)
curr.pnas.norm = normalizeVariance(curr.pnas.norm, details = TRUE)
curr.pnas.norm = log10(curr.pnas.norm$matnorm + 1)
curr.pnas.pca = RSpectra::svds(A = t(as.matrix(curr.pnas.norm)), k = 50,
                               opts = list(center = TRUE, scale = FALSE,
                                           maxitr = 2000, tol = 1e-10))
curr.pnas.pca = curr.pnas.pca$u
rownames(curr.pnas.pca) = rownames(pcs)

veloviz.pnas = buildVeloviz(
  curr = curr.pnas,
  proj = proj.pnas,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 50,
  similarity.threshold = 0.5,
  distance.weight = 0.01,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

#UMAP
par(mfrow = c(1,1))
set.seed(0)
nnGraph.pnas <- as_nn_graph(graph = veloviz.pnas$graph, k = 50)

emb.pnas.umap <- uwot::umap(X = NULL, nn_method = nnGraph.pnas, min_dist = 0.3)
rownames(emb.pnas.umap) = rownames(curr.pnas.pca)
plotEmbedding(emb.pnas.umap, colors = col, main = 'Xia et al cell cycle genes - UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")

