## Script for the basis processing of mIHC
##Yulan Deng, last updated 2025-12-26

###################################
#gating function for mIHC(R 4.3.1)#
###################################
markerAllPosGate_fun <- function(marker.label="CD45+CD3+CD8+",ROI_df=S2961ROI01_df,
marker.match=marker_match_df,cell.th=cellTh,ROI_stat=S2961ROI01_stat)
{
	marker.list <- unlist(strsplit(marker.label,"+",fixed=T))
	ROI_cell_df <- ROI_df[ROI_df[,"AreaShape_Area"]>=cell.th,]
	cb.label<- sapply(seq(length(marker.list)),function(x){
		res.cb.label <- paste(marker.list[seq(x)],collapse="+")
		return(res.cb.label)
	})
	cb.label <- paste0(cb.label,"+")
	No.cb.label <- ROI_stat[cb.label]
	ROI_tmp_df <- ROI_cell_df
	gate_th <- cell.th
	for(i in seq(length(marker.list)))
	{
		markerV <- ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]
		markerVsort <- sort(markerV,decreasing=T)
		if(No.cb.label[i]==0)
		{
			gate_tmp <- max(markerVsort)+1
		}else{
			gate_tmp <- markerVsort[No.cb.label[i]]
		}
		gate_th <- c(gate_th,gate_tmp)
		ROI_tmp_df <- ROI_tmp_df[ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]>=gate_tmp,]
	}
	names(gate_th) <- c("cell",marker.list)
	return(gate_th)	
}

Gate2No_fun <- function(marker.label="CD45+CD3+CD8-",gate.th=gate_test1,ROI_df=S2961ROI01_df,
marker.match=marker_match_df)
{
	marker.list1 <- unlist(strsplit(marker.label,"+",fixed=T))
	marker.list <- unlist(strsplit(marker.list1,"-",fixed=T))
	marker_index <- sapply(seq(length(marker.list)),function(x){
		indexTmp <- sum(nchar(marker.list[seq(x)]))+x
		indexRes <- substr(marker.label,indexTmp,indexTmp)
		return(indexRes)
	})
	
	ROI_cell_df <- ROI_df[ROI_df[,"AreaShape_Area"]>=gate.th[1],]
	ROI_tmp_df <- ROI_cell_df

	for(i in seq(length(marker.list)))
	{
		markerV <- ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]
		if(marker_index[i]=="+")
		{
			markerVsort <- sort(markerV,decreasing=T)
			ROI_tmp_df <- ROI_tmp_df[ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]>=gate.th[i+1],]
		}else{
			markerVsort <- sort(markerV,decreasing=F)
			ROI_tmp_df <- ROI_tmp_df[ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]<gate.th[i+1],]
		}
	}
	return(nrow(ROI_tmp_df))	
}

Gate2df_fun <- function(marker.label="CD45+CD3+CD8-",gate.th=Gate2No_test3,ROI_df=S2961ROI01_df,
marker.match=marker_match_df)
{
	marker.list1 <- unlist(strsplit(marker.label,"+",fixed=T))
	marker.list <- unlist(strsplit(marker.list1,"-",fixed=T))
	marker_index <- sapply(seq(length(marker.list)),function(x){
		indexTmp <- sum(nchar(marker.list[seq(x)]))+x
		indexRes <- substr(marker.label,indexTmp,indexTmp)
		return(indexRes)
	})
	
	ROI_cell_df <- ROI_df[ROI_df[,"AreaShape_Area"]>=gate.th[1],]
	ROI_tmp_df <- ROI_cell_df

	for(i in seq(length(marker.list)))
	{
		markerV <- ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]
		if(marker_index[i]=="+")
		{
			markerVsort <- sort(markerV,decreasing=T)
			ROI_tmp_df <- ROI_tmp_df[ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]>=gate.th[i+1],]
		}else{
			markerVsort <- sort(markerV,decreasing=F)
			ROI_tmp_df <- ROI_tmp_df[ROI_tmp_df[,marker.match[marker.match[,"Short"]==marker.list[i],"Full"]]<gate.th[i+1],]
		}
	}
	return(ROI_tmp_df)
}

subdf_pos_fun <- function(old.cb="CD45+CD3-",marker.list="CD20",ROI_df=CD45p_CD3n_sub_ROI_df,
marker.match=marker_match_df,ROI_stat=S2961ROI01_stat)
{
	cb.label <- paste0(old.cb,marker.list,"+")
	No.cb.label <- ROI_stat[cb.label]
	markerV <- ROI_df[,marker.match[marker.match[,"Short"]==marker.list,"Full"]]
	markerVsort <- sort(markerV,decreasing=T)
	if(No.cb.label==0)
	{
		gate_tmp <- max(markerVsort)+1
	}else{
		gate_tmp <- markerVsort[No.cb.label]
	}
	return(gate_tmp)	
}

subdf_2markers_fun <- function(old.cb="CD45+CD3+CD8+",marker1="PD1",marker2="EOMES",
ROI_df=CD45p_CD3p_CD8p_sub_ROI_df,
marker.match=marker_match_df,ROI_stat=S2961ROI01_stat)
{
	#marker1
	cb.label1 <- paste0(old.cb,marker1,"+",marker2,c("+","-"))
	No.cb.label1 <- sum(ROI_stat[cb.label1])
	markerV1 <- ROI_df[,marker.match[marker.match[,"Short"]==marker1,"Full"]]
	markerV1sort <- sort(markerV1,decreasing=T)
	if(No.cb.label1==0)
	{
		gate_tmp1 <- max(markerV1sort)+1
	}else{
		gate_tmp1 <- markerV1sort[No.cb.label1]
	}
	#marker2
	cb.label2 <- paste0(old.cb,marker1,c("+","-"),marker2,"+")
	No.cb.label2 <- sum(ROI_stat[cb.label2])
	markerV2 <- ROI_df[,marker.match[marker.match[,"Short"]==marker2,"Full"]]
	markerV2sort <- sort(markerV2,decreasing=T)
	if(No.cb.label2==0)
	{
		gate_tmp2 <- max(markerV2sort)+1
	}else{
		gate_tmp2 <- markerV2sort[No.cb.label2]
	}
	gate_tmp <- c(gate_tmp1,gate_tmp2)
	names(gate_tmp) <- c(marker1,marker2)
	return(gate_tmp)	
}
