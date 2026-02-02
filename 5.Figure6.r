## Script for the basis processing of visium HD
##Yulan Deng, last updated 2025-12-26

############################################
#spatial interaction for visium HD(R 4.4.0)#
############################################
library(RColorBrewer)
library(anndata)
library(Seurat)
library(dplyr)
library(patchwork)
library('harmony')
library(spatstat)
library(CellChat)
library(patchwork)

v2s_obj <- readRDS(file="v2s_obj.rds")
v2s_tumor <- readRDS(file="v2s_tumor.rds")
#spatial interaction
ad9 <- read_h5ad("953229-A_b2c.v4.h5ad")
ad8 <- read_h5ad("837317-3A_b2c.v4.h5ad")

cell_pos <- rbind(ad9$obs[,c("array_row","array_col")]*2,
ad8$obs[,c("array_row","array_col")]*2)
rownames(cell_pos) <- c(paste0(rownames(ad9$obs),"_1"),
paste0(rownames(ad8$obs),"_2"))
cellType <- c(as.character(ad9$obs[,"predicted_labels"]),
as.character(ad8$obs[,"predicted_labels"]))
names(cellType) <- rownames(cell_pos)
interCell <- intersect(rownames(cell_pos),rownames(v2s_tumor@meta.data))
cellType[interCell] <- paste0("Tumor_c",as.character(v2s_tumor@meta.data[,"RNA_snn_res.0.1"]))

pos_df <- data.frame(cellID=names(cellType),array_row=cell_pos[,1],array_col=cell_pos[,2],
predicted_labels=cellType,stringsAsFactors=F)

cell_cell_spatial_interaction_fun2 <- function(cell1="Tumor_c1",cell2="Stroma",roi_pos=pos_df)
{
	library(spatstat)
	roi.pos=roi_pos
	density_stat=table(roi_pos[,"predicted_labels"])
	if((density_stat[cell1]==0)|(density_stat[cell2]==0)){
		return(0)
	}else if(cell1!=cell2){
		cell1_index <- which(roi.pos[,"predicted_labels"]==cell1)
		cell2_index <- which(roi.pos[,"predicted_labels"]==cell2)
		cell1_obj <- ppp(x = roi.pos[cell1_index,"array_row"] ,
				   y = roi.pos[cell1_index,"array_col"],
				   window = owin(c(min(roi.pos[,"array_row"]), 
								   max(roi.pos[,"array_row"])),
								 c(min(roi.pos[,"array_col"]),
								   max(roi.pos[,"array_col"]))),
				   unitname = "micrometer")
		
		cell2_obj <- ppp(x = roi.pos[cell2_index,"array_row"] ,
				   y = roi.pos[cell2_index,"array_col"],
				   window = owin(c(min(roi.pos[,"array_row"]), 
								   max(roi.pos[,"array_row"])),
								 c(min(roi.pos[,"array_col"]),
								   max(roi.pos[,"array_col"]))),
				   unitname = "micrometer")
				   
		k1dist <- nncross(cell1_obj,cell2_obj,k=1,what = "dist")
		k1distMin <- min(k1dist)
		if(k1distMin>20)
		{
			return(0)
		}else{
			nk=min(length(cell2_index),50)
			kndist <- nncross(cell1_obj,cell2_obj,k=1:nk,what = "dist")
			return(sum(kndist<20)/(density_stat[cell1]*density_stat[cell2]))
		}
	}else if(cell1==cell2){
		cell1_index <- which(roi.pos[,"predicted_labels"]==cell1)
		cell2_index <- which(roi.pos[,"predicted_labels"]==cell2)
		if(length(cell1_index)==1)
		{
			return(0)
		}else{
			cell1_obj <- ppp(x = roi.pos[cell1_index,"array_row"] ,
					   y = roi.pos[cell1_index,"array_col"],
					   window = owin(c(min(roi.pos[,"array_row"]), 
									   max(roi.pos[,"array_row"])),
									 c(min(roi.pos[,"array_col"]),
									   max(roi.pos[,"array_col"]))),
					   unitname = "micrometer")
			
			cell2_obj <- ppp(x = roi.pos[cell2_index,"array_row"] ,
					   y = roi.pos[cell2_index,"array_col"],
					   window = owin(c(min(roi.pos[,"array_row"]), 
									   max(roi.pos[,"array_row"])),
									 c(min(roi.pos[,"array_col"]),
									   max(roi.pos[,"array_col"]))),
					   unitname = "micrometer")
				   
			k2dist <- nncross(cell1_obj,cell2_obj,k=2,what = "dist")
			k2distMin <- min(k2dist)
			if(k2distMin>20)
			{
				return(0)
			}else{
				nk=min(length(cell2_index),50)
				kndist <- nncross(cell1_obj,cell2_obj,k=2:nk,what = "dist")
				return((sum(kndist<20)/density_stat[cell1])/density_stat[cell2])
			}
		}
	}
}

cell_list <- setdiff(unique(cellType),paste0("Tumor_c",0:1))
cell_paird <- data.frame(cell1=rep(paste0("Tumor_c",0:1),length(cell_list)),
cell2=rep(cell_list,rep(2,length(cell_list))))

cell_cell_spatial_interaction_mt <- matrix(0,nrow=length(cell_list),ncol=2)
rownames(cell_cell_spatial_interaction_mt)=cell_list
colnames(cell_cell_spatial_interaction_mt)=paste0("Tumor_c",0:1)
for(i in seq(nrow(cell_paird)))
{
	res_tmp <- cell_cell_spatial_interaction_fun2(cell1=cell_paird[i,"cell1"],
	cell2=cell_paird[i,"cell2"],
	roi_pos=pos_df)
	cell_cell_spatial_interaction_mt[cell_paird[i,"cell2"],cell_paird[i,"cell1"]] <- res_tmp
	print(i)
}

save(cell_cell_spatial_interaction_mt,file="cell_cell_spatial_interaction_mt.visiumHD.RData")