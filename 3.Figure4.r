## Script for the basis processing of CyCIF
##Yulan Deng, last updated 2025-12-26

##################################
#tumor scale PCA by group kmeans)#
##################################
#smp为样本名
library(factoextra)
library(cluster)
library(RColorBrewer)
library(pheatmap)

set.seed(1)

#input
load(paste0(smp,".tumor_mt.2.RData"))

#scale
tumor_mt_scale <- scale(tumor_mt)

#PCA by group
comCS <- prcomp(tumor_mt[,c("Ecad", "H3K27me3" ,"Vimentin","p75" ,"MITF" )], center = TRUE,scale. = TRUE)
summCS<-summary(comCS)
nPCCS <- which(summCS$"importance"["Cumulative Proportion",]>=0.8)[1]

comD <- prcomp(tumor_mt[,c("pH2AX", "53BP1" )], center = TRUE,scale. = TRUE)

comP <- prcomp(tumor_mt[,c("CCNA2","PCNA","CCND1", "Ki67","pHH3" )], center = TRUE,scale. = TRUE)

comS <- prcomp(tumor_mt[,c("Axl", "PDGFR","pS6", "pGSK3b")], center = TRUE,scale. = TRUE)
summS<-summary(comS)
nPCS <- which(summS$"importance"["Cumulative Proportion",]>=0.8)[1]

f4 <- fviz_nbclust(cbind(comCS$"x"[,1:nPCCS],cbind(comD$"x"[,1],cbind(comP$"x"[,1],comS$"x"[,1:nPCS]))), kmeans, method = "silhouette", nstart = 10,iter.max = 20L)

Nk4 <- which(f4$"data"[,2]==max(f4$"data"[,2]))[1]

if(nrow(comCS$"x")<10000)
{
	ClassRes4 <- kmeans(cbind(comCS$"x"[,1:nPCCS],cbind(comD$"x"[,1],cbind(comP$"x"[,1],comS$"x"[,1:nPCS]))),centers=Nk4,iter.max = 20L,nstart = 10)
}else{
	sampleID <- sample(1:nrow(comCS$"x"),10000)
	ClassRes4 <- kmeans(cbind(comCS$"x"[sampleID,1:nPCCS],cbind(comD$"x"[sampleID,1],cbind(comP$"x"[sampleID,1],comS$"x"[sampleID,1:nPCS]))),centers=Nk4,iter.max = 20L,nstart = 10)
}
kmeansClass4 <- ClassRes4$"cluster"

save(kmeansClass4,file=paste0(smp,"_tumor_kmeans_best_sil.RData"))