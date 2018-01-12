# SNPs-based-Phylogeny
#The aim is to creat Maximum likelihood phylogenetic tree from a multiple alignment sequence file by using RStudio package of "phangorn" 

library(optrees)
library(vegan)
library(poppr)
library(pegas)
library(igraph)
install.packages("SyNet")
library(SyNet)
library(seedy)
MDR_ceppi.fasta <- system.file("/Users/ghodousi.arash/Desktop/NONE_joint_cf4_cr4_fr75_ph4_samples11_amended_u95_phylo_w12.plainIDs.fasta",package="adegenet")
obj.MDR_ceppi <- fasta2genlight("/Users/ghodousi.arash/Desktop/NONE_joint_cf4_cr4_fr75_ph4_samples11_amended_u95_phylo_w12.plainIDs.fasta")
MDR_ceppi.fasta <- fasta2DNAbin("/Users/ghodousi.arash/Desktop/NONE_joint_cf4_cr4_fr75_ph4_samples11_amended_u95_phylo_w12.plainIDs.fasta")
MDR_ceppi.dis <- dist.gene(MDR_ceppi.fasta, method = "pairwise",
                            variance = FALSE)
MDR_ceppi.matrix <- as.matrix(MDR_ceppi.dis, Class="graphNEL")
MDR_ceppi.Graph <- graph.adjacency(MDR_ceppi.matrix, mode ="undirected", weighted = TRUE )
MDR_ceppi.msn <- poppr.msn(obj.MDR_ceppi, MDR_ceppi.matrix,mlg.compute = "original",clustering.algorithm =  "nearest_neighbor", showplot = FALSE)
edgeWeights =E(MDR_ceppi.msn$graph)$weight
library("igraph")
set.seed(1234567890)
MDR_ceppi.msn <- plot_poppr_msn(obj.MDR_ceppi, MDR_ceppi.msn,gadj = 100,inds="All",gweight = 0,
                                 gscale=FALSE,vertex.label.cex=0.3,
                                 vertex.label.dist=0, nodebase = 1.25, layfun = layout.auto,
                                 pop.leg = TRUE,wscale=FALSE,edge.label=edgeWeights,
                                 edge.label.cex=0.3,edge.label.color="red",scale.leg = FALSE,
                                 edge.color="lightgrey",edge.width=0.5,
                                 edge.connectivity(MDR_ceppi.msn, source = NULL, target = NULL,checks=TRUE),
                                 vertex.color="white",
                                 vertex.shape="circle",vertex.frame.color="grey")
