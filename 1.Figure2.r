## Script for the basis processing of WES and nanostring
##Yulan Deng, last updated 2025-12-26

###############
#mutect (bash)#
###############
#${ref} is the file directory for reference genome
#${Tumor} is the name of tumor sample
#${Normal} is the name of normal sample
#${transcriptomeDir} is the directory for reference genome(hg38)
~/gatk/gatk-4.0.2.1/gatk Mutect2 -R $ref/Homo_sapiens_assembly19.fasta \
-I ${Tumor}.bam \
-I ${Normal}.bam \
--tumor-sample ${Tumor} -normal ${Normal} \
--TMP_DIR . --germline-resource $ref/af-only-gnomad.raw.sites.b37.vcf.gz \
--output ${Tumor}.vcf

#########################
#complexhetamap(R 4.0.3)#
#########################
library(gplots)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
load("loh_gene_df.RData")
load("mut_df.RData")
load("MHCclassI_genelist.RData")
load("cna_pt98_df.RData")
load("cna_pt42_pt208_df.RData")
load("MHCclassI_class_genelist.RData")

cna_df <- rbind(cna_pt42_pt208_df,cna_pt98_df)

#check sample, sample_df
sample_df <- data.frame(Usample=c(
"Pt208.9.10.14","Pt208.10.22.14","Pt208C1.3.11.15","Pt208C2.3.11.15",
"Pt208.5.13.15","Pt208.8.19.15","Pt42.10.17.14","Pt42.11.24.14",
"Pt42.12.24.14","Pt42.12.29.14","Pt98.1.25.17","Pt98.2.5.14", 
"Pt98.3.3.14","Pt98.9.7.16"),
SampleID1=c("SRR6504374","SRR6504375","SRR6504376","SRR6504377",
"SRR6504378","SRR6504379","SRR6504398","SRR6504399",
"SRR6504400","SRR6504401","Pt98.1.25.17","Pt98.2.5.14", 
"Pt98.3.3.14","Pt98.9.7.16"),
SampleID2=c("SRR6504374","SRR6504375","SRR6504376","SRR6504377",
"SRR6504378","SRR6504379","SRR6504398","SRR6504399",
"SRR6504400","SRR6504401","Pt98_1.25.17","Pt98_2.5.14", 
"Pt98_3.3.14","Pt98.9.7.16"),stringsAsFactors=F)

setdiff(cna_df[,"Sample"],sample_df[,2])
setdiff(loh_gene_df[,"Sample"],sample_df[,2])
setdiff(mut_df[,"Sample"],sample_df[,3])
sample_order <- c(
"Pt208.9.10.14","Pt208.10.22.14","Pt208C1.3.11.15","Pt208C2.3.11.15",
"Pt208.5.13.15","Pt208.8.19.15","Pt42.10.17.14","Pt42.11.24.14",
"Pt42.12.24.14","Pt42.12.29.14","Pt98.2.5.14", 
"Pt98.3.3.14","Pt98.9.7.16","Pt98.1.25.17")

gene_sel <- intersect(unique(c(mut_df[,"gene"],loh_gene_df[,"Gene.Name"])),MHCclassI_genelist)

gene_stat <- sort(table(c(mut_df[mut_df[,"gene"]%in%gene_sel,"gene"],
loh_gene_df[loh_gene_df[,"Gene.Name"]%in%gene_sel,"Gene.Name"])),decreasing=T)
geneF <- names(gene_stat)[gene_stat>2]

geneF_list <- lapply(MHCclassI_class_genelist,function(x) intersect(x,geneF))
sapply(geneF_list,length)
geneF1 <-  geneF_list["Antigen processing-Cross presentation"][[1]]
geneF2 <- setdiff(geneF_list["Antigen Presentation: Folding, assembly and peptide loading of class I MHC"][[1]],geneF1)
geneF3 <- setdiff(geneF_list["Antigen processing: Ubiquitination & Proteasome degradation"][[1]],c(geneF1,geneF2))
geneF_order <- c(geneF1,geneF2,geneF3)
col_ha <- HeatmapAnnotation(subPathway=c(rep("c1",length(geneF1)),
rep("c2",length(geneF2)),
rep("c3",length(geneF3))),col=list(subPathway=c("c1"=
brewer.pal(5, 'Set2')[1],"c2"=brewer.pal(5, 'Set2')[2],"c3"=brewer.pal(5, 'Set2')[3])))

discrete_mat <- matrix(0,nrow=length(sample_order),ncol=length(geneF_order))
colnames(discrete_mat) <- geneF_order
rownames(discrete_mat) <- sample_order

discrete_mat1 <-discrete_mat
mutM_df <- mut_df[mut_df[,"gene"]%in%geneF_order,]
mutM_df[,"Sample"] <- sample_df[match(mutM_df[,"Sample"],sample_df[,3]),1]
mutM_df[mutM_df[,"ExonicFunc.refGene"]==".","ExonicFunc.refGene"] <- "SpliceSite"
mutM_df[mutM_df[,"ExonicFunc.refGene"]=="nonsynonymous SNV","ExonicFunc.refGene"] <- "Missense"
mutM_df[mutM_df[,"ExonicFunc.refGene"]=="stopgain","ExonicFunc.refGene"] <- "Nonsense"
mutM_df[mutM_df[,"ExonicFunc.refGene"]=="frameshift substitution","ExonicFunc.refGene"] <- "FrameShiftSub"
for(i in seq(nrow(mutM_df)))
{
	discrete_mat1[mutM_df[i,"Sample"],mutM_df[i,"gene"]] <- mutM_df[i,"ExonicFunc.refGene"]
}

discrete_mat2 <-discrete_mat
cnaM_df <- cna_df[cna_df[,"gene"]%in%geneF_order,]
cnaM_df[,"Sample"] <- sample_df[match(cnaM_df[,"Sample"],sample_df[,2]),1]
for(i in seq(nrow(cnaM_df)))
{
	discrete_mat2[cnaM_df[i,"Sample"],cnaM_df[i,"gene"]] <- cnaM_df[i,"status"]
}

discrete_mat3 <-discrete_mat
lohM_df <- loh_gene_df[loh_gene_df[,"Gene.Name"]%in%geneF_order,]
lohM_df[,"Sample"] <- sample_df[match(lohM_df[,"Sample"],sample_df[,2]),1]
for(i in seq(nrow(lohM_df)))
{
	discrete_mat3[lohM_df[i,"Sample"],lohM_df[i,"Gene.Name"]] <- "LossOfHeterozygosity"
}

discrete_colors <- c("0" = "#FFFFFF", "Nonsense" = "#F18080", "Missense" = "#447DBF", 
                    "SpliceSite" = "#90C73E", "FrameShiftSub" = "#FFE000")
names(discrete_colors) <- c("0","Nonsense","Missense","SpliceSite","FrameShiftSub")

stat_snp <- t(apply(discrete_mat1,1,function(x) {
	res_1 <- sum(x=="Nonsense")
	res_2 <- sum(x=="Missense")
	res_3 <- sum(x=="SpliceSite")
	res_4 <- sum(x=="FrameShiftSub")
	res <- c(res_1,res_2,res_3,res_4)
	return(res)
}))
stat_cna <- t(apply(discrete_mat2,1,function(x) {
	res_1 <- sum(x=="Deletion")
	res_2 <- sum(x=="Amplification")
	res <- c(res_1,res_2)
	return(res)
}))
stat_loh <- apply(discrete_mat3,1,function(x) sum(x=="LossOfHeterozygosity"))

row_ha <- rowAnnotation(
  "NucleotideVariants" = anno_barplot(stat_snp, 
                       gp = gpar(fill = c("#F18080","#447DBF","#90C73E","#FFE000")),
                       height = unit(2, "cm")),
  "CopyNumberVariants" = anno_barplot(stat_cna, 
                     gp = gpar(fill = NA,col = c("#7075B7","#F37575")),
                     height = unit(2, "cm")),
  "LOH" = anno_barplot(stat_loh, 
                         gp = gpar(fill = NA,col="#000000"),
                         height = unit(2, "cm")),
  gap = unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10)
)

border_col <- matrix("#DDDDDD", nrow(discrete_mat2), ncol(discrete_mat2))
border_col[discrete_mat2=="Deletion"] <- "#7075B7"
border_col[discrete_mat2=="Amplification"] <- "#F37575"

pdf(file="genomic_alteration_of_antigen_presentation_pathway_rotation_group.pdf",
width=16,height=5)
Heatmap(discrete_mat1,
        name = "Nucleotide variations", 
        col = discrete_colors, 
        row_gap = unit(0, "mm"),
        column_gap = unit(0, "mm"),
		
		cluster_rows=FALSE,
		cluster_columns=FALSE,

        rect_gp = gpar(type = "none"),
		
		right_annotation=row_ha,
		top_annotation=col_ha,
        
        layer_fun = function(j, i, x, y, width, height, fill) {

            gap <- unit(1.5, "mm")
            
            cell_width <- width - gap
            cell_height <- height - gap
			
            grid.rect(x = x, y = y, 
                     width = cell_width, 
                     height = cell_height,
                     gp = gpar(fill = fill, col = border_col,lwd=2))

			show_text <- discrete_mat3=="LossOfHeterozygosity"
			grid.text(ifelse(show_text,sprintf("%s","*"),""), x, y, gp = gpar(fontsize = 20))
			
        })
dev.off()

######################################
#nanostring Pt35 RNA cluster(R 4.3.1)#
######################################
library(readr)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(cluster)

nanostring_RNA <- read_tsv(file="Geomean normalized_Gao_RNA samples only_all data.txt")
Pt35_group <- read_tsv(file="Pt35_Group.txt")

nanostring_RNA_df <- as.data.frame(nanostring_RNA)
Pt35_group_df <- as.data.frame(Pt35_group)

##check names
all(Pt35_group_df[,1]%in%colnames(nanostring_RNA_df))

nanostring_mt <- t(as.matrix(nanostring_RNA_df[,-1]))
colnames(nanostring_mt) <- nanostring_RNA_df[,1]
submt <- nanostring_mt[Pt35_group_df[,1],]

submt_log2 <- log2(submt+1)
rownames(submt_log2) <- Pt35_group_df[,"Note_for_biopsies"]

submt_log2_scale <- t(scale(submt_log2))
#check row gene,col samples
submt_log2_scale[1:5,1:5]
f4 <- fviz_nbclust(submt_log2_scale, kmeans, method = "silhouette",  nstart = 10,iter.max = 20L)

pdf(file="PT35.nanostring.RNA.cluster.best.k.pdf")
f4
dev.off()

Nk4 <- which(f4$"data"[,2]==max(f4$"data"[,2]))[1]


ClassRes4 <- kmeans(submt_log2_scale,centers=Nk4)
kmeansClass4 <- ClassRes4$"cluster"
save(ClassRes4 ,file="ClassRes4.PT35.nanotringRNA.RData")
kmeansClass4_f <- factor(as.character(kmeansClass4),levels=as.character(sort(unique(kmeansClass4))))
submt_log2_order <- submt_log2[,order(kmeansClass4_f)]

annotation_col = data.frame(Cluster = sort(kmeansClass4_f))
rownames(annotation_col) = colnames(submt_log2_order)

pt_ht <- pheatmap(submt_log2_order, 
               color = colorRampPalette(c('green','black','red'))(100), 
               border_color = NA,  
               scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
			   annotation_col = annotation_col,
               show_rownames = TRUE, 
               show_colnames = FALSE)

pdf("PT35.nanotringRNA.cluster.pdf",width=14,height=3)
pt_ht
dev.off()

