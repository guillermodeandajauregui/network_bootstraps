library("igraph")
library("data.table")

readSIF <- function(SIFfile){
  #read a SIF file into R 
  g<-fread(SIFfile, data.table = FALSE)
  g<-g[order(g$V2, decreasing = TRUE),] #ensure ordering
  names(g)[names(g)=="V2"] <- "weight"
  g<-g[,c(1,3,2)]
  return(g)
}

genesymclean <- function(g){
  #remove links that include non-gene symbol probes from Affymetrix chips 
  g_clean <- g[!grepl(pattern = "AFF*", x = g$V1) & !grepl(pattern = "AFF*", x = g$V3),]
  g_clean <- g_clean[!grepl(pattern = "*_at", x = g_clean$V1) & !grepl(pattern = "*_at", x = g_clean$V3),]
  return(g_clean)
}

sif2cleannx <- function(SIFfile){
  return(graph.data.frame(d = genesymclean(readSIF(SIFfile), directed = FALSE)))
}

pvalue<-function(mi, n=100){
  alfa = 1.062
  beta = -48.7
  gamma = -0.634
  p = exp(alfa -mi*(-beta + (-gamma * n)))
  return(p)
}

pdf.sif<-function(sif){
  return(density(x = sif$weight))
}

plot_density.wo<-function(..., file = "plot.pdf"){
  dots = list(...)
  colores = rainbow(n = length(seq_along(dots)))
  pdf(file = file, paper = "a4r")
  plot(x = dots[[1]]$x,
       y = dots[[1]]$y,
       xlab = "MI",
       ylab = "Frequency",
       col = colores[1],
       type = "l")
  
  for(i in seq_along(dots)){
    if(i > 1){
      points(x = dots[[i]]$x,
             y = dots[[i]]$y,
             type = "l",
             col = colores[i])
    }
  }
  dev.off()
}

nodesfromcompsize <- function(clusters, size){
  gps <- groups(clusters)
  nclist <- list()
  for (i in gps){
    if(length(i) >= size){
      nclist <-c(nclist, i)}
  }
  return(unlist(nclist))
}

nw_analysis<-function(nw){
  nodes_nwc <- length(V(nw))
  edges_nwc <- length(E(nw))
  clusters_nwc <- clusters(nw)$no
  clusteringcoefficient_nwc <- transitivity(nw)
  max_weight <- max(E(nw)$weight)
  min_weight <- min(E(nw)$weight)
  oput<-rbind(nodes_nwc, edges_nwc, clusters_nwc, clusteringcoefficient_nwc, max_weight, min_weight)
  rownames(oput)<-c("nodes", "edges", "clusters", "cc", "max_weight", "min_weight")
  return(oput)
}

nw_analysis_p<-function(nw, N=100){
  nodes_nwc <- length(V(nw))
  edges_nwc <- length(E(nw))
  clusters_nwc <- clusters(nw)$no
  clusteringcoefficient_nwc <- transitivity(nw)
  max_weight <- max(E(nw)$weight)
  min_weight <- min(E(nw)$weight)
  pval <- pvalue(mi = min_weight, n = N)
  oput<-rbind(nodes_nwc, edges_nwc, clusters_nwc, clusteringcoefficient_nwc, max_weight, min_weight, pval)
  rownames(oput)<-c("nodes", "edges", "clusters", "cc", "max_weight", "min_weight", "pvalue")
  return(oput)
}