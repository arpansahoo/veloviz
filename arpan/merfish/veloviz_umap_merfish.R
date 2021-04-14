graphics.off()

library(veloviz)
library(velocyto.R)
library(uwot)

#### LOAD MERFISH DATA ####
#load from veloviz package
col <- MERFISH$col
pcs <- MERFISH$pcs
vel <- MERFISH$vel

curr <- vel$current
proj <- vel$projected

merfish.genes = rownames(curr)

nuc <- MERFISH$nuc
cyto <- MERFISH$cyto
all <- nuc + cyto

#### MERFISH CELL CYCLE GENES TO CALCULATE PSEUDOTIME ####
#MERFISH genes exhibiting cell-cycle-dependent expression (Xia et al 2019, Supp Dataset 8)
# https://www.pnas.org/content/116/39/19490
cycle.genes.pnas <- read.csv("pnas.01.csv",header = TRUE)
cycle.genes.pnas <- cycle.genes.pnas$Gene

merfish.cycle.pnas <- intersect(merfish.genes, cycle.genes.pnas)
nuc.pnas <- nuc[merfish.cycle.pnas,]
cyto.pnas <- cyto[merfish.cycle.pnas,]
all.pnas <- nuc.pnas + cyto.pnas

curr.pnas <- curr[merfish.cycle.pnas,]
proj.pnas <- proj[merfish.cycle.pnas,]


#veloviz
veloviz.pnas <-  buildVeloviz(
  curr = curr.pnas,
  proj = proj.pnas,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 50,
  similarity.threshold = -1,
  distance.weight = 0.01,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

par(mfrow = c(1,1))
emb.pnas.vv <- veloviz.pnas$fdg_coords

n.cells <- nrow(emb.pnas.vv)

emb.vv.scaled <- scale(emb.pnas.vv, scale=F)
emb.vv.scaled[,1] <- emb.vv.scaled[,1]/max(abs(emb.vv.scaled[,1]))
emb.vv.scaled[,2] <- emb.vv.scaled[,2]/max(abs(emb.vv.scaled[,2]))
par(mfrow = c(1,1))

ptime.vv <- 180*(atan(emb.vv.scaled[,2]/emb.vv.scaled[,1]))/pi 
names(ptime.vv) = rownames(emb.pnas.vv)

neg.x <- which(emb.vv.scaled[,1]<0)
ptime.vv[neg.x] <- -1*ptime.vv[neg.x]
ptime.vv <- (ptime.vv + 90)/180 

col.gradient = colorRampPalette(c('blue','red'))
ptime.col <- col.gradient(n.cells)[as.numeric(cut(ptime.vv,breaks = n.cells))]
names(ptime.col) = rownames(emb.pnas.vv)


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
par(mfrow=c(1,1), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))
# plotEmbedding(emb.all.vv, col = ptime.col[rownames(emb.all.vv)],
#               main = 'VeloViz', xlab = "VeloViz X", ylab = "VeloViz Y")
# 
# # Normal UMAP
# set.seed(0)
# emb.all.umap = uwot::umap(X = pcs[,1:3], min_dist = 0.3)
# rownames(emb.all.umap) = rownames(pcs)
# plotEmbedding(emb.all.umap, col = ptime.col,
#               main = 'UMAP', xlab = "UMAP X", ylab = "UMAP Y")

# UMAP (initialized with velo)
set.seed(0)
nnGraph.all <- veloviz::asNNGraph(veloviz.all)
emb.all.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph.all, min_dist = 0.3)
rownames(emb.all.umapVelo) <- rownames(emb.all.vv)
plotEmbedding(emb.all.umapVelo, col = ptime.col[rownames(emb.all.umapVelo)],
              main = 'UMAP with Veloviz',
              xlab = "UMAP X", ylab = "UMAP Y")


# show velocities
pdf("merfish_new.pdf")

# par(mfrow=c(2,2), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))
# show.velocity.on.embedding.cor(emb.all.vv, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.all.vv)],0.4),
#                                cex = 1, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 2,
#                                grid.n = 30, arrow.lwd = 1.25, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
#                                main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Y")
# 
# show.velocity.on.embedding.cor(emb.all.umap, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col,0.4),
#                                cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 2,
#                                grid.n = 30, arrow.lwd = 1.25, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
#                                main = "UMAP", xlab = "UMAP X", ylab = "UMAP Y")

show.velocity.on.embedding.cor(emb.all.umapVelo, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.all.umapVelo)],0.4),
                               cex = 1, arrow.scale = 0.75, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "UMAP with VeloViz", xlab = "UMAP X", ylab = "UMAP Y")

dev.off()