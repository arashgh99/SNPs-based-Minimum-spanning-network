install.packages("poppr")
library(poppr)


My.fasta <- system.file("/Path-to-fasta",package="adegenet")
My.obj <- fasta2genlight("/Path-to-fasta") # to creat genlight objects which store SNP data efficiently by packing binary allele calls into single bits
My.fasta <- fasta2DNAbin("/Path-to-fasta")
My.dis <- dist.gene(My.fasta, method = "pairwise",
                            variance = FALSE).    #To compute a matrix of distances between pairs of individuals from a matrix or a data frame of genetic data
My.matrix <- as.matrix(My.dis, Class="graphNEL")
My.Graph <- graph.adjacency(My.matrix, mode ="undirected", weighted = TRUE ) #for creating igraph graphs from adjacency matrices
My.msn <- poppr.msn(My.obj, My.matrix,mlg.compute = "original",clustering.algorithm =  "nearest_neighbor", showplot = FALSE) # to Create minimum spanning network of selected populations

library("igraph")
set.seed()
My.msn <- plot_poppr_msn(My.obj, My.msn,gadj = 100,inds="All",gweight = 0,
                                 gscale=FALSE,vertex.label.cex=0.3,
                                 vertex.label.dist=0, nodebase = 1.25, layfun = layout.auto,
                                 pop.leg = TRUE,wscale=FALSE,
                                 edge.label.cex=0.3,edge.label.color="red",scale.leg = FALSE,
                                 edge.color="lightgrey",edge.width=0.5,
                                 edge.connectivity(MDR_ceppi.msn, source = NULL, target = NULL,checks=TRUE),
                                 vertex.color="white",
                                 vertex.shape="circle",vertex.frame.color="grey") #To plot the minimum spanning network