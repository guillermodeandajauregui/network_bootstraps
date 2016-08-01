# Bootstrapped validation of networks
library("igraph")
library("data.table")
#Step 1: What is the network structure
#a) For each subtype, plot whole Mutual Information distribution.

  #read nw
g<-fread("BASAL_20000.sif", data.table = FALSE)
g<-g[order(g$V2, decreasing = TRUE),]
#b) Identify min MI value for the 10000 interaction
min_mi<- min(g$V2)
#c) Calculate equivalent p-value for number of samples and that p-value
pvalue<-function(mi, n){
alfa = 1.062
beta = -48.7
gamma = -0.634
p = exp(alfa -mi*(-beta + (-gamma * n)))
return(p)
}

min_pvalue = pvalue(min_mi)
#igraph calculate values
nw<-graph_from_data_frame(d = g, directed = FALSE)

nodes <- length(V(nw))
clusters <- clusters(nw)
clusteringcoefficient<- transitivity(nw, typ)

##Remove non-GS rows 
g_clean <- g[!grepl(pattern = "AFF*", x = g$V1) & !grepl(pattern = "AFF*", x = g$V3),]
g_clean <- g_clean[!grepl(pattern = "*_at", x = g_clean$V1) & !grepl(pattern = "*_at", x = g_clean$V3),]


#Step 2: Repeat for bootstrapped networks
rows = length((g[[1]]))
iters = c(rows*.01, rows*.1, rows*.9)


# top k, top 100k, top million.
#when single component is found, go through that interval until percolation. 
#maybe try the "cleaning less than 20 sized islands" thing before
#maybe networkx for this job?