## Script for the basis processing of DSP
##Yulan Deng, last updated 2025-12-26

#######################################################
#Calculation of MHCII-enriched region in DSP (R 4.4.1)#
#######################################################
library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(ggpubr)
DSP_raw <- read_tsv(file="DSP20240627T.txt")
DSP_df <- as.data.frame(DSP_raw)
DSP_mt <- as.matrix(DSP_df[,-1])
rownames(DSP_mt) <- DSP_df[,1]
sample_label <- read_tsv(file="DSP20250507_label.txt")
sample_label_df <- as.data.frame(sample_label)
sample_label_NAomit_df <- na.omit(sample_label_df)
DSP_mt_scale <- scale(DSP_mt)
protein_list <- c("CD4","CD8","CD3","HLA-DR","CD11c")
MHCII_score <- apply(DSP_mt_scale[,protein_list],1,mean)
#plot
#NR pre on  post
smp_list1 <- sample_label_NAomit_df[(sample_label_NAomit_df[,"CLASS"]=="NR")&(sample_label_NAomit_df[,"Status"]=="Pre"),"NAME"]
smp_list2 <- sample_label_NAomit_df[(sample_label_NAomit_df[,"CLASS"]=="NR")&(sample_label_NAomit_df[,"Status"]=="On"),"NAME"]
smp_list3 <- sample_label_NAomit_df[(sample_label_NAomit_df[,"CLASS"]=="NR")&(sample_label_NAomit_df[,"Status"]=="Post"),"NAME"]

res_df <- data.frame(MHCIIscore=c(MHCII_score[smp_list1],MHCII_score[smp_list2],
MHCII_score[smp_list3]),
TreatmentTime=factor(c(rep("Pre",length(smp_list1)),
rep("On",length(smp_list2)),
rep("Post",length(smp_list3))),levels=c("Pre","On","Post")))

my_comparisons <- list()
my_comparisons[[1]] <- c("Post","On")
my_comparisons[[2]] <- c("Pre","Post")

pdf(file="DSP_MHCII_score.NR.pre.on.post.20251203.pdf")
p1 <- ggplot(res_df, aes(x = TreatmentTime, y = MHCIIscore)) +
  geom_boxplot(outlier.shape = NA, aes(color = TreatmentTime)) +  # 绘制箱线图
  geom_jitter(width = 0.2, aes(color = TreatmentTime)) +  # 添加抖动散点
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
theme(panel.background=element_rect(fill='transparent', color='black'),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
  text=element_text(color ="black"),
  axis.text.x=element_text(color ="black"),
  axis.text.y=element_text(color ="black"))+ # 
  labs(title = "MHCII score in NR(DSP)",
       x = "Treatment time",
       y = "MHCII score")+scale_fill_manual(values =brewer.pal(8, "Set1")[3:1] ) 
print(p1)

dev.off()

#post NR vs R
my_comparisons <- list()
my_comparisons[[1]] <- c("NR","R")

smp_list1 <- sample_label_NAomit_df[(sample_label_NAomit_df[,"CLASS"]=="NR")&(sample_label_NAomit_df[,"Status"]=="Post"),"NAME"]
smp_list2 <- sample_label_NAomit_df[(sample_label_NAomit_df[,"CLASS"]=="R")&(sample_label_NAomit_df[,"Status"]=="Post"),"NAME"]

res_df <- data.frame(MHCIIscore=c(MHCII_score[smp_list1],MHCII_score[smp_list2]),
Response=factor(c(rep("NR",length(smp_list1)),
rep("R",length(smp_list2)))))

pdf(file="DSP_MHCIIscore_score.post.NRvsR.20251203.pdf")
p1 <- ggplot(res_df, aes(x = Response, y = MHCIIscore)) +
  geom_boxplot(outlier.shape = NA, aes(color = Response)) +  # 绘制箱线图
  geom_jitter(width = 0.2, aes(color = Response)) +  # 添加抖动散点
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
theme(panel.background=element_rect(fill='transparent', color='black'),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
  text=element_text(color ="black"),
  axis.text.x=element_text(color ="black"),
  axis.text.y=element_text(color ="black"))+ # 
  labs(title = "MHCII score in post-PD1 (DSP)",
       x = "Response",
       y = "MHCII score")+scale_fill_manual(values =brewer.pal(8, "Set2") ) 
print(p1)

dev.off()
