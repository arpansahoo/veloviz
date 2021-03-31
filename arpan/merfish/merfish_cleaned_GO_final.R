#### SETUP ####
library(veloviz)
library(velocyto.R)
library(uwot)
source("../as_nn_graph.R")

#### LOAD MERFISH DATA ####
#load from veloviz package
col <- MERFISH$col
# pcs <- MERFISH$pcs
vel <- MERFISH$vel

curr <- vel$current
proj <- vel$projected

merfish.genes = rownames(curr)

nuc <- MERFISH$nuc
cyto <- MERFISH$cyto
all <- nuc + cyto

#### MERFISH CELL CYCLE GENES FOR PSEUDOTIME ####
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
plotEmbedding(emb.pnas.vv, colors = col[rownames(emb.pnas.vv)],
              main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Y")

n.cells <- nrow(emb.pnas.vv)

emb.vv.scaled <- scale(emb.pnas.vv, scale=F)
emb.vv.scaled[,1] <- emb.vv.scaled[,1]/max(abs(emb.vv.scaled[,1]))
emb.vv.scaled[,2] <- emb.vv.scaled[,2]/max(abs(emb.vv.scaled[,2]))
par(mfrow = c(1,1))
plot(emb.vv.scaled, pch=16)

ptime.vv <- 180*(atan(emb.vv.scaled[,2]/emb.vv.scaled[,1]))/pi 
names(ptime.vv) = rownames(emb.pnas.vv)

neg.x <- which(emb.vv.scaled[,1]<0)
ptime.vv[neg.x] <- -1*ptime.vv[neg.x]
ptime.vv <- (ptime.vv + 90)/180 

# ptime.vv[which(ptime.vv<0)] <- -1*ptime.vv[which(ptime.vv<0)]


col.gradient = colorRampPalette(c('blue','red'))
ptime.col <- col.gradient(n.cells)[as.numeric(cut(ptime.vv,breaks = n.cells))]
names(ptime.col) = rownames(emb.pnas.vv)

plot(scale(emb.pnas.vv), col=ptime.col, pch=16)


#### GO CELL CYCLE GENES (GO:0000278) ####
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_MITOTIC_CELL_CYCLE - no longer accessible? renamed to below
# cycle.genes.go = read.csv("GO_0000278.csv",header = FALSE)
# cycle.genes.go = cycle.genes.go$V1
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_MITOTIC_CELL_CYCLE.html
cycle.genes.go <- read.csv("geneset.csv", header = FALSE)
cycle.genes.go <- as.matrix(cycle.genes.go)

merfish.cycle.go <- intersect(merfish.genes, cycle.genes.go)
nuc.go <- nuc[merfish.cycle.go,]
cyto.go <- cyto[merfish.cycle.go,]
all.go <- nuc.go + cyto.go

#normalize and reduce dims 
all.go.norm = normalizeDepth(all.go)
all.go.norm = normalizeVariance(all.go.norm)
all.go.norm = log10(all.go.norm + 1)
pcs.go <- reduceDimensions(all.go.norm, center = T, scale = T)

curr.go <- curr[merfish.cycle.go,]
proj.go <- proj[merfish.cycle.go,]

veloviz.go = buildVeloviz(
  curr = curr.go,
  proj = proj.go,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 65, #60, #45, 
  similarity.threshold = -1, 
  distance.weight = 0.1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.go.vv = veloviz.go$fdg_coords

par(mfrow = c(2,2))
#VV
plotEmbedding(emb.go.vv, colors = ptime.col[rownames(emb.go.vv)],
              main = 'GO cell cycle genes - veloviz')

#UMAP
set.seed(0)
emb.go.umap = umap::umap(pcs.go[,1:3], min_dist = 0.3)$layout
rownames(emb.go.umap) = rownames(pcs.go)
plotEmbedding(emb.go.umap, colors = ptime.col, main = 'GO cell cycle genes - UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")


pdf("merfish_full_GOccGenes.pdf")
par(mfrow=c(2,2), omi = c(0.1,0.1,0.1,0.1))
show.velocity.on.embedding.cor(emb.go.vv, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.vv)],0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Y",
                               cex.main = 2, cex.lab = 1.5)
show.velocity.on.embedding.cor(emb.go.umap, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.umap)],0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = 'UMAP', xlab = "UMAP X", ylab = "UMAP Y",
                               cex.main = 2, cex.lab = 1.5)
dev.off()

#### G2M SCORE CELLS ####
#get cell cycle genes
cc.genes <- read.csv("cell_cycle_genes.csv")

#convert gene names
# https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
library(AnnotationHub)
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
library(dplyr)
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
library(GenomicFeatures)
library(ensembldb)
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")
# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cc.genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")


# # https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION.html
# g2m_genes <- read.csv("./full_data/GO_0044839_g2m.txt", header = FALSE)
# g2m_genes <- as.matrix(g2m_genes)

# normalize all counts 
all.norm <- normalizeDepth(all)
all.norm <- normalizeVariance(all.norm)
all.norm <- log10(all.norm + 1)
all.norm = t(scale(t(as.matrix(all.norm))))
all.norm[all.norm < 0] = 0

# G2M genes in normalized MERFISH data 
g2m.merfish.genes <- intersect(rownames(all.norm), g2m_genes)

#plot marker gene expression
par(mfrow=c(4,4), mar=c(1,1,1,1))
for (g in g2m.merfish.genes){
  print(g)
  gene.exp = all.norm[g,]
  # hist(gene.exp, breaks=300, main=g)
  col.gradient = colorRampPalette(c('white','red'))
  gene.col = col.gradient(n.cells)[as.numeric(cut(gene.exp,breaks = n.cells))]
  names(gene.col) = rownames(emb.pnas.vv)
  
  plotEmbedding(emb.pnas.vv, colors = gene.col, main=g, xlab="", ylab="")
}

#G2M score = sum of individual genes 
g2m.sum <- colSums(all.norm[g2m.merfish.genes,])
g2m.sum <- g2m.sum/max(g2m.sum)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
hist(g2m.sum, breaks=300, main="g2m gene sum")
col.gradient = colorRampPalette(c('white','red'))
gene.col = col.gradient(n.cells)[as.numeric(cut(g2m.sum,breaks = n.cells))]
names(gene.col) = rownames(emb.pnas.vv)
plotEmbedding(emb.pnas.vv, colors = gene.col, xlab="", ylab="", alpha = 1, cex = 1)


#remove G2M cells 
# pdf("merfish_g2m_score.pdf", width = 10)
par(mfrow=c(1,2))
hist(g2m.sum, breaks=600, main = "Cell G2M Scores", xlab = "G2M Score", ylab = "Count", cex.main = 2,cex.lab = 1.5)
abline(v=0.3, lwd = 2, col = "red")
plot(emb.pnas.vv, col = gene.col, pch=16,
     xlab = "VeloViz X", ylab = "VeloViz Y", main = "G2M Cells", xaxt = "n", yaxt="n",
     cex.main = 2, cex.lab = 1.5)

#remove cells with g2m score
not.g2m.cells <- colnames(curr)[-which(g2m.sum>0.3)]
points(emb.pnas.vv[-which(rownames(emb.pnas.vv) %in% not.g2m.cells),], pch = 4, col = "blue", cex = 0.7, lwd = 0.5)
legend(x=5, y=-3, legend = c("G2M score","G2M cell"), pch = c(16,4), col = c("red", "blue"), pt.cex = c(1,0.7))
# dev.off()

#remove 90% of g2m cells
g2m.cells <- colnames(curr)[which(g2m.sum>0.3)]
set.seed(1)
g2m.cells.keep <- sample(g2m.cells,length(g2m.cells)/10)
not.g2m.cells <- c(not.g2m.cells, g2m.cells.keep)

par(mfrow = c(1,2))
plot(emb.go.vv, col = gene.col, pch=16,
     xlab = "VeloViz X", ylab = "VeloViz Y", main = "G2M Cells", xaxt = "n", yaxt="n",
     cex.main = 2, cex.lab = 1.5)
plot(emb.go.vv[not.g2m.cells,], col = gene.col[not.g2m.cells], pch=16,
     xlab = "VeloViz X", ylab = "VeloViz Y", main = "G2M Cells", xaxt = "n", yaxt="n",
     cex.main = 2, cex.lab = 1.5)

#### EMBEDDINGS W/O SUBSET - GO GENES ####
#recalculate velocity with subset of data
cells.sub <- not.g2m.cells
nuc.sub <- nuc[,cells.sub]
cyto.sub <- cyto[,cells.sub]

#normalize
all.sub <- nuc.sub + cyto.sub

cpm <- normalizeDepth(all.sub)
varnorm <- normalizeVariance(cpm)
lognorm <- log10(varnorm + 1)

pcs.sub <- reduceDimensions(lognorm, center = T, scale=T, nPCs=50)
cell.dist <- as.dist(1-cor(t(pcs.sub)))

vel.sub = gene.relative.velocity.estimates(cyto.sub,
                                           nuc.sub,
                                           kCells = 30,
                                           cell.dist = cell.dist,
                                           fit.quantile = 0.1)
curr.sub <- vel.sub$current
proj.sub <- vel.sub$projected

merfish.gene.sub <- rownames(curr.sub)
merfish.cycle.go.sub <- intersect(merfish.gene.sub, cycle.genes.go)

curr.go.sub <- curr.sub[merfish.cycle.go.sub,]
proj.go.sub <- proj.sub[merfish.cycle.go.sub,]

veloviz.go.sub = buildVeloviz(
  curr = curr.go.sub,
  proj = proj.go.sub,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  pca = TRUE,
  nPCs = 3,
  center = TRUE,
  scale = TRUE,
  k = 45, 
  similarity.threshold = 0, 
  distance.weight = 0.05,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.go.vv = veloviz.go.sub$fdg_coords

par(mfrow = c(2,2))
#VV
plotEmbedding(emb.go.vv, colors = ptime.col[rownames(emb.go.vv)],
              main = 'GO cell cycle genes - veloviz')

#PCA
emb.go.pca = pcs.sub[,1:2]
plotEmbedding(emb.go.pca, colors = ptime.col[rownames(emb.go.pca)], main = 'GO cell cycle genes - pca')

#tSNE
set.seed(0)
emb.go.tsne = Rtsne::Rtsne(pcs.sub[,1:3], perplexity = 50)$Y
rownames(emb.go.tsne) = rownames(pcs.sub)
plotEmbedding(emb.go.tsne, colors = ptime.col[rownames(emb.go.tsne)], main = 'GO cell cycle genes - t-SNE',
              xlab = "t-SNE X", ylab = "t-SNE y")

#UMAP
set.seed(0)
emb.go.umap = umap::umap(pcs.sub[,1:3], min_dist = 0.3)$layout
rownames(emb.go.umap) = rownames(pcs.sub)
plotEmbedding(emb.go.umap, colors = ptime.col[rownames(emb.go.umap)], main = 'GO cell cycle genes - UMAP',
              xlab = "UMAP X", ylab = "UMAP Y")


pdf("merfish_gap_GOccGenes.pdf")
par(mfrow=c(2,2), omi = c(0.1,0.1,0.1,0.1), mai = c(0.82,0.82,0.62,0.22))
show.velocity.on.embedding.cor(emb.go.vv, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.vv)],0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Y",
                               cex.main = 2, cex.lab = 1.5)
show.velocity.on.embedding.cor(emb.go.pca, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.pca)],0.4),
                               cex = 1, arrow.scale = 10, show.grid.flow = TRUE, min.grid.cell.mass = 1, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = 'PCA', xlab = "PC 1", ylab = "PC 2",
                               cex.main = 2, cex.lab = 1.5)
show.velocity.on.embedding.cor(emb.go.tsne, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.tsne)],0.4),
                               cex = 1, arrow.scale = 2, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = 't-SNE', xlab = "t-SNE X", ylab = "t-SNE Y",
                               cex.main = 2, cex.lab = 1.5)
show.velocity.on.embedding.cor(emb.go.umap, vel, n = 100, scale = 'sqrt', cell.colors = scales::alpha(ptime.col[rownames(emb.go.umap)],0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 2, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = 'UMAP', xlab = "UMAP X", ylab = "UMAP Y",
                               cex.main = 2, cex.lab = 1.5)
dev.off()
