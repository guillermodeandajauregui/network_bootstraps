source("bootstrapping_functions.R")
#load all interactions
sif <- genesymclean(readSIF(SIFfile = "test_case.sif"))
#plot MI distribution
plot_density.wo(pdf.sif(sif), file = "mi_density.pdf")
#generate network
nw<-graph_from_data_frame(d = sif, directed = FALSE)
#analyze network
analysis<-nw_analysis(nw)
write.csv(analysis, file = "analysis_total.txt")

#repeat for 90%, 10%, 1%, 0.1%, 0.01%, 0.001%
interactionstotal = length((sif[[1]]))
intot<-interactionstotal
iters = c(interactionstotal*.9, interactionstotal*.1, interactionstotal*.01, interactionstotal*.001, interactionstotal*.0001)
rm(nw) #for better memory usage

##90%
sif90 <- sif[1:round(intot*.90),]
nw90 <- graph_from_data_frame(d = sif90, directed = FALSE)
analysis<-nw_analysis(nw90)
write.csv(analysis, file = "analysis_90.txt")
###90%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw90), size = 20)
nw90ns <- induced_subgraph(nw90, V(nw90)[V(nw90)$name%in%ncomps])
analysis<-nw_analysis(nw90ns)
write.csv(analysis, file = "analysis_90ns.txt")
rm(nw90)
rm(nw90ns)

##10%
sif10 <- sif[1:round(intot*.10),]
nw10 <- graph_from_data_frame(d = sif10, directed = FALSE)
analysis<-nw_analysis(nw10)
write.csv(analysis, file = "analysis_10.txt")
###10%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw10), size = 20)
nw10ns <- induced_subgraph(nw10, V(nw10)[V(nw10)$name%in%ncomps])
analysis<-nw_analysis(nw10ns)
write.csv(analysis, file = "analysis_10ns.txt")
rm(nw10)
rm(nw10ns)

##01%
sif01 <- sif[1:round(intot*.01),]
nw01 <- graph_from_data_frame(d = sif01, directed = FALSE)
analysis<-nw_analysis(nw01)
write.csv(analysis, file = "analysis_01.txt")
###01%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw01), size = 20)
nw01ns <- induced_subgraph(nw01, V(nw01)[V(nw01)$name%in%ncomps])
analysis<-nw_analysis(nw01ns)
write.csv(analysis, file = "analysis_01ns.txt")
rm(nw01)
rm(nw01ns)

##001%
sif001 <- sif[1:round(intot*.001),]
nw001 <- graph_from_data_frame(d = sif001, directed = FALSE)
analysis<-nw_analysis(nw001)
write.csv(analysis, file = "analysis_001.txt")
###001%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw001), size = 20)
nw001ns <- induced_subgraph(nw001, V(nw001)[V(nw001)$name%in%ncomps])
analysis<-nw_analysis(nw001ns)
write.csv(analysis, file = "analysis_001ns.txt")
rm(nw001)
rm(nw001ns)

##0001%
sif0001 <- sif[1:round(intot*.0001),]
nw0001 <- graph_from_data_frame(d = sif0001, directed = FALSE)
analysis<-nw_analysis(nw0001)
write.csv(analysis, file = "analysis_0001.txt")
###0001%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw0001), size = 20)
nw0001ns <- induced_subgraph(nw0001, V(nw0001)[V(nw0001)$name%in%ncomps])
analysis<-nw_analysis(nw0001ns)
write.csv(analysis, file = "analysis_0001ns.txt")
rm(nw0001)
rm(nw0001ns)

##00001%
sif00001 <- sif[1:round(intot*.00001),]
nw00001 <- graph_from_data_frame(d = sif00001, directed = FALSE)
analysis<-nw_analysis(nw00001)
write.csv(analysis, file = "analysis_00001.txt")
###00001%, no small islands (ns)
ncomps <- nodesfromcompsize(clusters(nw00001), size = 20)
nw00001ns <- induced_subgraph(nw00001, V(nw00001)[V(nw00001)$name%in%ncomps])
analysis<-nw_analysis(nw00001ns)
write.csv(analysis, file = "analysis_00001ns.txt")
rm(nw00001)
rm(nw00001ns)