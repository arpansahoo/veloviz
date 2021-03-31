#### SETUP AND BUILD GRAPH ####
library(veloviz)
library(velocyto.R)
library(igraph)
library(plot3D)
library(viridis)


clusters = pancreas$clusters # cell type annotations
pcs = pancreas$pcs # PCs used to make other embeddings (UMAP, tSNE..)
vel = pancreas$vel # velocity

#choose colors based on clusters for plotting later
# cell.cols = rainbow(8)[as.numeric(clusters)]
# names(cell.cols) = names(clusters)
col = rev(plasma(length(levels(clusters))))
cell.cols = col[clusters] 
names(cell.cols) = names(clusters)

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

emb.veloviz = veloviz$fdg_coords
par(mfrow = c(1,1))
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='veloviz')


g <- plotVeloviz(veloviz, clusters=clusters)


#### PLOT 3D FR ####
dim3 <- igraph::layout_with_fr(veloviz$graph, dim=3)
rownames(dim3) <- colnames(curr)

pdf("panc_3D_1.pdf", width = 8, height = 4)
par(mfrow=c(1,3))
plot(dim3[,1:2], pch=16, col=cell.cols, xlab="VeloViz X", ylab="VeloViz Y", xaxt = "n", yaxt="n")
plot(dim3[,c(1,3)], pch=16, col=cell.cols, xlab="VeloViz X", ylab="VeloViz Z", xaxt = "n", yaxt="n")
plot(dim3[,c(2,3)], pch=16, col=cell.cols, xlab="VeloViz Y", ylab="VeloViz Z", xaxt = "n", yaxt="n")
dev.off()


colvar <- as.numeric(clusters)
col.3d <- unique(cell.cols)
col.3d[unique(colvar)] <- col.3d


par(mfrow=c(1,3))
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=0, phi=0, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=0, phi=90, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=90, phi=0, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)

pdf("panc_3D_2.pdf", width = 8, height = 4)
par(mfrow=c(1,3))
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=45, phi=0, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=70, phi=0, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=-20, phi=10, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)
dev.off()

par(mfrow = c(1,1))
plot3D::points3D(dim3[,1], dim3[,2], dim3[,3], theta=70, phi=20, pch=16, cex=0.7, colvar=colvar, col=col.3d, colkey=FALSE)


## plot arrows 
par(mfrow = c(1,3))
show.velocity.on.embedding.cor(scale(dim3[,1:2]), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Y")
show.velocity.on.embedding.cor(scale(dim3[,2:3]), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "VeloViz", xlab = "VeloViz X", ylab = "VeloViz Z")
show.velocity.on.embedding.cor(scale(dim3[,c(2,3)]), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "VeloViz", xlab = "VeloViz Y", ylab = "VeloViz Z")

#### OTHER LAYOUTS ####

#Davidson and Harel - nope
dh <- igraph::layout_with_dh(veloviz$graph, coords = pcs[,1:2])
rownames(dh) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(dh, groups=clusters[rownames(dh)], main='veloviz - dh')



#GEM - nope
gem <- igraph::layout_with_gem(veloviz$graph, coords = pcs[,1:2])
rownames(gem) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(gem, groups=clusters[rownames(gem)], main='veloviz - gem')


#Graphopt --good
graphopt <- igraph::layout_with_graphopt(veloviz$graph)
rownames(graphopt) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(graphopt, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - graphopt')



#Kamada Kawai - nope
kk <- igraph::layout_with_kk(veloviz$graph)
rownames(kk) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(kk, groups=clusters[rownames(kk)], main='veloviz - kk')



#Large graph layout -- good
lgl <- igraph::layout_with_lgl(veloviz$graph)
rownames(lgl) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(lgl, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - lgl')


#Multidimensional layout - nope 
mds <- igraph::layout_with_mds(veloviz$graph)
rownames(mds) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(mds, groups=clusters[rownames(mds)], main='veloviz - mds')



#Sugiyama layout - nope 
sugiyama <- igraph::layout_with_sugiyama(veloviz$graph)
rownames(sugiyama$layout) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(sugiyama$layout, groups=clusters[rownames(sugiyama$layout)], main='veloviz - sugiyama')



#DRL layout
drl <- igraph::layout_with_drl(veloviz$graph)
rownames(drl) <- colnames(curr)

par(mfrow = c(1,2))
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(drl, colors = cell.cols[rownames(emb.veloviz)], main='veloviz - drl')


## good layouts 

par(mfrow = c(2,2))
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz - fr')
plotEmbedding(graphopt, groups=clusters[rownames(graphopt)], main='veloviz - graphopt')
plotEmbedding(lgl, groups=clusters[rownames(lgl)], main='veloviz - lgl')
plotEmbedding(drl, groups=clusters[rownames(drl)], main='veloviz - drl')


par(mfrow = c(2,2))
show.velocity.on.embedding.cor(scale(emb.veloviz), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "veloviz - fr", xlab = "VeloViz X", ylab = "VeloViz Y")
show.velocity.on.embedding.cor(scale(graphopt), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "veloviz - graphopt", xlab = "VeloViz X", ylab = "VeloViz Y")
show.velocity.on.embedding.cor(scale(lgl), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "veloviz - lgl", xlab = "VeloViz X", ylab = "VeloViz Y")
show.velocity.on.embedding.cor(scale(drl), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1.5, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "veloviz - drl", xlab = "VeloViz X", ylab = "VeloViz Y")


pdf("panc_layouts.pdf", width = 8, height = 4)
par(mfrow = c(1,3))
plotEmbedding(emb.veloviz, colors = cell.cols[rownames(emb.veloviz)], main='Fruchterman and Reingold')
plotEmbedding(graphopt, colors = cell.cols[rownames(emb.veloviz)], main='GraphOpt')
#https://arxiv.org/abs/2007.03619
plotEmbedding(lgl, colors = cell.cols[rownames(emb.veloviz)], main='Large Graph Layout')
#http://www.marcottelab.org/paper-pdfs/jmb-lgl.pdf
dev.off()


pdf("panc_layouts_vel.pdf", width = 8, height = 4)
par(mfrow = c(1,1))
show.velocity.on.embedding.cor(scale(emb.veloviz), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "Fruchterman and Reingold", xlab = "", ylab = "")
show.velocity.on.embedding.cor(scale(graphopt), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "GraphOpt", xlab = "", ylab = "")
show.velocity.on.embedding.cor(scale(lgl), vel, n = 50, scale = 'sqrt', cell.colors = scales::alpha(cell.cols,0.4),
                               cex = 1, arrow.scale = 1, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, 
                               grid.n = 30, arrow.lwd = 1, do.par = FALSE, frame.plot = TRUE, xaxt = 'n', yaxt = 'n',
                               main = "Large Graph Layout", xlab = "", ylab = "")
dev.off()

