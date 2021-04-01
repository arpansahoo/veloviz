#### SETUP ####
library(velocyto.R)
library(veloviz)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
source("as_nn_graph.R")
load("10X_mouse_dimred.RData")

#### RUN WITH SUBSAMPLED DATA ####

sample.size <- c(100, 500, 1000, 2500, 5000, 7500, 9295) #runtimes
sampled.cells <- list()

for (n in sample.size){
  curr.size.samples <- list()

  for (r in c(1:3)){
    set.seed(r)
    curr.cells <- sample(cells.keep, n)
    curr.size.samples[[r]] <- curr.cells
  }

  sampled.cells[[as.character(n)]] <- curr.size.samples
}


vv.times <- matrix(0, nrow = length(sample.size), ncol = 3)
umap.times <- matrix(0, nrow = length(sample.size), ncol = 3)
umapVelo.times <- matrix(0, nrow = length(sample.size), ncol = 3)

for (n in c(1:length(sample.size))){
  print(paste("n =", n))
  curr.size.samples <- sampled.cells[[n]]
  

  for (r in c(1:3)){
    print(paste("r = ",r))
    curr.cells <- curr.size.samples[[r]]
    cell.dist <- as.dist(1-cor(t(pcs[curr.cells,])))

    curr.s <- s[,curr.cells]
    curr.u <- u[,curr.cells]

    vel <- gene.relative.velocity.estimates(curr.s,curr.u,kCells = 30, cell.dist = cell.dist, fit.quantile = 0.1)

    curr <- vel$current
    proj <- vel$projected

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
      similarity.threshold = 0.1, 
      distance.weight = 1,
      distance.threshold = 0.6,
      weighted = TRUE,
      seed = 0,
      verbose = FALSE
    )

    vv.start <- Sys.time()
    g <- plotVeloviz(veloviz, clusters=com, seed = 0)
    vv.end <- Sys.time()

    vv.time <- difftime(vv.end, vv.start, units = "secs")
    vv.times[n,r] <- vv.time

    ##UMAP
    set.seed(0)
    emb.umap = uwot::umap(pcs[curr.cells,], min_dist = 0.5, ret_nn = TRUE)

    umap.start <- Sys.time()
    set.seed(0)
    emb.umapBackIn = uwot::umap(X = NULL, nn_method = emb.umap$nn, min_dist = 0.5)
    umap.end <- Sys.time()

    umap.time <- difftime(umap.end, umap.start, units = "secs")
    umap.times[n,r] <- umap.time
    
    ##UMAP with VeloViz

    # Get NN graph
    nnGraph <- as_nn_graph(graph = veloviz$graph, k = 15)

    # input nnGraph to UMAP and plot
    set.seed(0)
    umapVelo.start <- Sys.time()
    emb.umapVelo <- uwot::umap(X = NULL, nn_method = nnGraph, min_dist = 0.5)
    umapVelo.end <- Sys.time()
    umapVelo.time <- difftime(umapVelo.end, umapVelo.start, units = "secs")
    umapVelo.times[n,r] <- umapVelo.time
  }
  
}

save(vv.times, umap.times, umapVelo.times, file = "runtimes.RData")

#### EVALUATE RUNTIME ####

sample.size <- c(100, 500, 1000, 2500, 5000, 7500, 9295)
load("runtimes.RData")

# VV RUNTIME #
rownames(vv.times) <- sample.size
vv.times <- t(vv.times)
vv.times.avg <- colMeans(vv.times)

vv.times.df <- data.frame(vv.times)
vv.times.df <- gather(vv.times.df)
vv.times.df$key <- unlist(lapply(sample.size, function(x) rep(x,3)))
colnames(vv.times.df) <- c("Number of Cells", "Runtime (seconds)")

vvt.avg.df <- data.frame(vv.times.avg)
vvt.avg.df <- cbind(sample.size, vvt.avg.df)
colnames(vvt.avg.df) <- c("Number of Cells", "Average Runtime (seconds)")


figtheme <-  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust=0.5, size = 36),
        axis.line=element_line(size=1),
        legend.position = 'none')

p <- ggplot(vv.times.df, aes(x=`Number of Cells`, y=`Runtime (seconds)`)) + 
  geom_point(colour = "red", size=4, shape=4) + 
  stat_summary(fun=mean, geom="line", colour="black", size = 1) +
  figtheme

p

ggsave(file = "vv_scalability.pdf", plot = p, width = 4, height = 4)


# UMAP RUNTIME #
rownames(umap.times) <- sample.size
umap.times <- t(umap.times)
umap.times.avg <- colMeans(umap.times)

umap.times.df <- data.frame(umap.times)
umap.times.df <- gather(umap.times.df)
umap.times.df$key <- unlist(lapply(sample.size, function(x) rep(x,3)))
colnames(umap.times.df) <- c("Number of Cells", "Runtime (seconds)")

umap.avg.df <- data.frame(umap.times.avg)
umap.avg.df <- cbind(sample.size, umap.avg.df)
colnames(umap.avg.df) <- c("Number of Cells", "Average Runtime (seconds)")


figtheme <-  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust=0.5, size = 36),
        axis.line=element_line(size=1),
        legend.position = 'none')

p <- ggplot(umap.times.df, aes(x=`Number of Cells`, y=`Runtime (seconds)`)) + 
  geom_point(colour = "red", size = 4, shape=4) + 
  stat_summary(fun=mean, geom="line", colour="black", size = 1) +
  figtheme

p

ggsave(file = "umap_scalability.pdf", plot = p, width = 4, height = 4)


# UMAP with Velo RUNTIME #
rownames(umapVelo.times) <- sample.size
umapVelo.times <- t(umapVelo.times)
umapVelo.times.avg <- colMeans(umapVelo.times)

umapVelo.times.df <- data.frame(umapVelo.times)
umapVelo.times.df <- gather(umapVelo.times.df)
umapVelo.times.df$key <- unlist(lapply(sample.size, function(x) rep(x,3)))
colnames(umapVelo.times.df) <- c("Number of Cells", "Runtime (seconds)")

umapVelo.avg.df <- data.frame(umapVelo.times.avg)
umapVelo.avg.df <- cbind(sample.size, umapVelo.avg.df)
colnames(umapVelo.avg.df) <- c("Number of Cells", "Average Runtime (seconds)")


figtheme <-  theme_bw(base_size = 18) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(hjust=0.5, size = 36),
        axis.line=element_line(size=1),
        legend.position = 'none')

p <- ggplot(umapVelo.times.df, aes(x=`Number of Cells`, y=`Runtime (seconds)`)) + 
  geom_point(colour = "red", size = 4, shape=4) + 
  stat_summary(fun=mean, geom="line", colour="black", size = 1) +
  figtheme

p

ggsave(file = "umap_velo_scalability.pdf", plot = p, width = 4, height = 4)

