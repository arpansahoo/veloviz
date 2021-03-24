library(velocyto.R)
library(veloviz)
library(hdf5r)
library(loomR)
library(enrichR)
library(MERINGUE)
library(plot3D)
library(viridis)
library(reticulate)
##load data
load("10X_mouse_unfiltered.RData")
load("10X_mouse_filtered.RData")
load("10X_mouse_dimred.RData")
col.clust <- rainbow(length(table(com)))[as.numeric(com)]
# col.clust <- plasma(length(table(com)))[as.numeric(com)]
names(com) <- names(col.clust) <- cells.keep

# lfile <- connect(filename = "SC3_v3_NextGem_SI_Neuron_10K_possorted_genome_bam_QLREO.loom", mode = "r+", skip.validate = T)
# counts <- lfile[["matrix"]]$read()
# genes <- lfile[["row_attrs"]][["Gene"]]$read()
# cells <- lfile[["col_attrs"]][["CellID"]]$read()
# cells <- sapply(cells, function(x) substr(x, 57, 72))
# 
# s <- lfile[["layers"]][["spliced"]]$read()
# u <- lfile[["layers"]][["unspliced"]]$read()
# # ambiguous <- lfile[["layers"]][["ambiguous"]]$read()
# # spanning <- lfile[["layers"]][["spanning"]]$read()
# 
# colnames(s) <- colnames(u) <- genes
# rownames(s) <- rownames(u) <- cells
# 
# s <- t(s)
# u <- t(u)
# all <- s+u
# 
# s.means <- rowMeans(s)
# u.means <- rowMeans(u)
# plot(log10(s.means+1),log10(u.means+1))
# 
# hist(log10(rowSums(all)+1), breaks=1000, ylim=c(0,200))
# hist(log10(colSums(all)+1), breaks=1000)
# 
# #save unfiltered
# save(s, u, cells, genes, file="10X_mouse_unfiltered.RData")
load("10X_mouse_unfiltered.RData")


#filter
# par(mfrow=c(2,2))
# hist(log10(rowSums(s)+1), breaks=1000, ylim=c(0,200))
# hist(log10(colSums(s)+1), breaks=1000)
# 
# hist(log10(rowSums(u)+1), breaks=1000, ylim=c(0,200))
# hist(log10(colSums(u)+1), breaks=1000)
# 
# genes.keep <- genes[which(log10(rowSums(s)+1)>2 & log10(rowSums(u)+1)>2)]
# cells.keep <- cells[which(log10(colSums(s)+1)>2 & log10(colSums(s)+1)>2)]
# 
# # genes.keep <- genes[which(log10(rowSums(all)+1)>2.5)]
# # cells.keep <- cells[which(log10(colSums(all)+1)>3.25)]
# 
# s <- s[genes.keep,cells.keep]
# u <- u[genes.keep,cells.keep]
# 
# save(s, u, cells.keep, genes.keep, file="10X_mouse_filtered.RData")
load("10X_mouse_filtered.RData")


# all <- all[genes.keep,cells.keep]
all <- s+u

s.means <- rowMeans(s)
u.means <- rowMeans(u)
par(mfrow=c(1,1))
plot(log10(s.means+1),log10(u.means+1))

hist(log10(rowSums(all)+1), breaks=1000)
hist(log10(colSums(all)+1), breaks=1000)

#there are duplicate rownames?? 
# genes.dup <- genes.keep[which(duplicated(genes.keep))]
# genes.keep <- genes.keep[which(!(genes.keep %in% genes.dup))]

# s <- s[genes.keep,]
# u <- u[genes.keep,]
# all <- all[genes.keep,]

#normalize
cpm <- normalizeDepth(all)
varnorm <- normalizeVariance(cpm)
lognorm <- log10(varnorm + 1)

#pca
pcs <- reduceDimensions(lognorm, center = T, scale = T, nPCs = 50)

par(mfrow=c(2,2))
emb.pca <- pcs[,1:2]
plot(emb.pca, pch=16)

#tSNE
set.seed(0)
emb.tsne = Rtsne::Rtsne(pcs, perplexity=30)$Y
rownames(emb.tsne) = rownames(pcs)
plot(emb.tsne, main='tSNE',pch = 16,
     xlab = "t-SNE X", ylab = "t-SNE Y")

##UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plot(emb.umap, main='UMAP',pch = 16,
     xlab = "UMAP X", ylab = "UMAP Y")

##diffusion map 
set.seed(0)
diffmap = destiny::DiffusionMap(pcs)
emb.diffmap = destiny::eigenvectors(diffmap)[,1:2]
row.names(emb.diffmap) = cells.keep
plot(emb.diffmap, pch = 16, main = "Diffusion map")

### clustering 
com <- MUDAN::getComMembership(pcs, 50, igraph::cluster_louvain)
col.clust <- rainbow(length(table(com)))[as.numeric(com)]
# col.clust <- plasma(length(table(com)))[as.numeric(com)]
names(com) <- names(col.clust) <- cells.keep

par(mfrow=c(2,2))
plot(emb.pca, pch=16, col=col.clust, cex=0.7, main="pca")
plot(emb.tsne, pch=16, col=col.clust, cex=0.7, main="tsne")
plot(emb.umap, pch=16, col=col.clust, cex=0.7, main="umap")
plot(emb.diffmap, pch=16, col=col.clust, cex=0.7, main="diffusion map")


#save for later
# save(s, u, lognorm, cells.keep, genes.keep, pcs, com, file="10X_mouse_dimred.RData")
load("10X_mouse_dimred.RData")


##velocity
cell.dist <- as.dist(1-cor(t(pcs)))

# vel <- gene.relative.velocity.estimates(s,u,kCells = 30, cell.dist = cell.dist, fit.quantile = 0.1)
# saveRDS(vel, file = "10_mouse_vel_01.rds")
vel <- readRDS(file = "10_mouse_vel_01.rds")
curr <- vel$current
proj <- vel$projected


vel.genes <- rownames(curr)
par(mfrow=c(1,1))
plot(log10(rowSums(s[vel.genes,]+1)), log10(rowSums(u[vel.genes,]+1)), pch = 16)


# #process curr and proj for veloviz 
# ##normalize depth 
# curr.norm = normalizeDepth(curr)
# proj.norm = normalizeDepth(proj)
# # #one cell gives NaNs?? 
# # curr.norm <- curr.norm[,-which(colnames(curr.norm) %in% c("AGGGCCTGTAGGACCA"))]
# # proj.norm <- proj.norm[,-which(colnames(curr.norm) %in% c("AGGGCCTGTAGGACCA"))]
# 
# 
# ##variance stabilize
# curr.varnorm.info = normalizeVariance(curr.norm, details = TRUE)
# curr.varnorm = curr.varnorm.info$matnorm
# 
# #use same model for projected
# scale.factor = curr.varnorm.info$df$scale_factor #gene scale factors
# names(scale.factor) = rownames(curr.varnorm.info$df)
# 
# m = proj.norm
# rmean = Matrix::rowMeans(m) #row mean
# sumx = Matrix::rowSums(m)
# sumxx = Matrix::rowSums(m^2)
# rsd = sqrt((sumxx - 2 * sumx * rmean + ncol(m) * rmean ^ 2) / (ncol(m)-1)) #row sd
# 
# proj.varnorm = proj.norm / rsd * scale.factor[names(rsd)]
# proj.varnorm = proj.norm[rownames(curr.varnorm),]



veloviz <- buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 50,
  center = TRUE,
  scale = TRUE,
  k = 30,
  similarity.threshold = 0,
  distance.weight = 1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.veloviz = veloviz$fdg_coords
par(mfrow=c(1,1))
plot(emb.veloviz, pch=16)


g <- plotVeloviz(veloviz, clusters=com, col=col.clust)


par(mfrow=c(1,1))
show.velocity.on.embedding.cor(emb.umap, vel, 
                               n = 100,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')
show.velocity.on.embedding.cor(emb.veloviz, vel, 
                               n = 100,
                               scale='rank',
                               cex=1, arrow.scale=2, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')


### annotation 
#get differential genes 
dg <- MUDAN::getDifferentialGenes(lognorm, com)
dg.sig <- lapply(dg, function(x) {
  #rownames(x)[x$Z > 8]
  x <- x[x$highest,]
  x <- x[x$Z > 4,]
  x <- x[order(x$Z, decreasing=TRUE),]
  na.omit(rownames(x))
})
length(unique(unlist(dg.sig)))

dg.sig


##gsea
#get mouse lists
dbs <- listEnrichrDbs()
dbs <- c("Mouse_Gene_Atlas")


#gsea 
enr <- list()

for (c in names(dg.sig)){
  sig.c <- enrichr(dg.sig[[c]], dbs)
  sig.c <- sig.c[[1]][1:10,1:4]
  enr[[c]] <- sig.c
}

enr

g = 'CDK1' ## Cycling
gexp <- scale(matnorm0[g,])[,1]
t = 2
gexp[gexp < -t] <- -t
gexp[gexp > t] <- t
plotEmbedding(emb.all[good.cells.scrnaseq,], col=MERINGUE:::map2col(gexp), main=g, alpha=0.1)


#### TESTING VELOVIZ ####
veloviz0 <- buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 50,
  center = TRUE,
  scale = TRUE,
  k = 30,
  similarity.threshold = 0,
  distance.weight = 1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.veloviz0 = veloviz0$fdg_coords
plot(emb.veloviz0, pch=16, col=col.clust, cex=0.7)

####*****
veloviz <- buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 50,
  center = TRUE,
  scale = TRUE,
  k = 50,
  similarity.threshold = 0,
  distance.weight = 1,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.veloviz = veloviz$fdg_coords
plot(emb.veloviz, pch=16, col=col.clust, cex=0.7)


emb.veloviz0 = veloviz0$fdg_coords
plot(emb.veloviz0, pch=16, col=col.clust, cex=0.7)

###### this one
veloviz <- buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 50,
  center = TRUE,
  scale = TRUE,
  k = 15,
  similarity.threshold = 0.1, #,
  distance.weight = 1,
  distance.threshold = 0.6, #1
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)


# emb.veloviz0 = veloviz0$fdg_coords
# plot(emb.veloviz, pch=16, col=col.clust, cex=0.7)

par(mfrow=c(1,1))
emb.veloviz = veloviz$fdg_coords
plot(emb.veloviz, pch=16, col=col.clust[rownames(emb.veloviz)], cex=0.7)


### use scvelo to visualize 
library(reticulate)
use_condaenv("cellrank", required=T)
scv <- import("scvelo")
plt  = import('matplotlib.pyplot')
sc = import("scanpy")

a <- scv$read("SC3_v3_NextGem_SI_Neuron_10K_possorted_genome_bam_QLREO.loom")
a$obs_names <- reticulate::r_to_py(cells)
a <- a[a$obs_names$isin(reticulate::r_to_py(cells.keep))]
a$obs['cluster'] <- reticulate::r_to_py(com)
# a2 <- a
# a2 <- a2[,a2$var_names$isin(reticulate::r_to_py(genes.keep))]

# a <- a[a$var_names$isin(reticulate::r_to_py(genes.keep))]
# a$var_names_make_unique()
# a$obs_names_make_unique()

scv$pp$filter_and_normalize(a, min_shared_counts = as.integer(20), n_top_genes = as.integer(2000))
sc$tl$pca(a)
sc$pp$neighbors(a, n_pcs = as.integer(30), n_neighbors = as.integer(30))
sc$tl$umap(a)
sc$tl$tsne(a)
sc$tl$diffmap(a)

scv$pp$moments(a, n_pcs = 30, n_neighbors = 30)

scv$tl$recover_dynamics(a)
scv$tl$velocity(a, mode = "dynamical")
scv$tl$velocity_graph(a)

# scv$tl$latent_time(a.vel)

##save for later
py_save_object(a, filename = "10Xmouse_embryo_scvelo_vel02.rds")
##read saved velocity
a = py_load_object("10Xmouse_embryo_scvelo_vel02.rds")

#plot embeddings 

sc$pl$pca(a, color='cluster')
sc$pl$umap(a, color='cluster')
sc$pl$tsne(a, color='cluster')
sc$pl$diffmap(a, color='cluster')




### 









#get colors 
plot.cols = unique(col.clust)
col = rev(rainbow(length(levels(com))))
c = unlist(sapply(c(1:length(col)), function(x) which(plot.cols == col[x])))
plot.cols = plot.cols[c]


#plot
pt.size = 50
dnsty = 0.9
plt.size = c(3.5,3.5)

## pca 
ax = scv$pl$velocity_embedding_stream(a, basis='umap',
                                      density = dnsty, cutoff_perc = 20, n_neighbors = 40L, title = "", size = pt.size, 
                                      show = FALSE, figsize = plt.size, palette = plot.cols, legend_loc = "on data")
# plt$show()
plt$savefig("full_umap.png")




















veloviz1 <- buildVeloviz(
  curr = curr, proj = proj,
  normalize.depth = TRUE,
  use.ods.genes = FALSE,
  alpha = 0.05,
  pca = TRUE,
  nPCs = 50,
  center = TRUE,
  scale = TRUE,
  k = 15,
  similarity.threshold = 0.2,
  distance.weight = 10,
  distance.threshold = 1,
  weighted = TRUE,
  seed = 0,
  verbose = FALSE
)

emb.veloviz1 = veloviz1$fdg_coords
plot(emb.veloviz1, pch=16, col=col.clust[row.names(emb.veloviz1)], cex=0.7)

### plot velocity 
show.velocity.on.embedding.cor(emb.umap, vel, 
                               n = 100,
                               scale='rank',
                               cex=1, arrow.scale=4, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')

show.velocity.on.embedding.cor(emb.veloviz, vel, 
                               n = 100,
                               scale='rank',
                               cex=1, arrow.scale=4, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')

show.velocity.on.embedding.cor(emb.veloviz, vel, 
                               n = 100,
                               scale='log',
                               cex=1, arrow.scale=4, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')

show.velocity.on.embedding.cor(emb.veloviz, vel, 
                               n = 100,
                               scale='sqrt',
                               cex=1, arrow.scale=4, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')

show.velocity.on.embedding.cor(emb.veloviz, vel, 
                               n = 100,
                               scale='linear',
                               cex=1, arrow.scale=4, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=1,do.par = F,
                               cell.colors=col.clust, main='VeloViz')



### plot 3d ###

plot(dim3[,1:2], pch=16, col=col.clust)
plot(dim3[,c(1,3)], pch=16, col=col.clust)
plot(dim3[,c(2,3)], pch=16, col=col.clust)

dim3 <- igraph::layout_with_fr(veloviz$graph, dim=3)
rownames(dim3) <- colnames(curr)
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=0, phi=0, pch=16, cex=0.7, colvar=as.numeric(com), col=(unique(col.clust)), colkey=FALSE)

                 