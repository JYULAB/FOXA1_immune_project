w18_DEG <- read.table(file = "~/R_project/FOXA1_luly_gain_loss/data/w18_DEG.csv",header = T, sep = ",")
w18_all_FOXA1 <- read.csv(file = "~/R_project/FOXA1_luly_gain_loss/data/FOXA1_peak_spike_in_xout_all_anno.csv",header = T)
w18_all_peakscore <- read.csv(file = "~/R_project/FOXA1_luly_gain_loss/data/FOXA1_peak_spike_in_xout_all.csv",header = T)
TSS_site <- read.csv("~/R_project/FOXA1_luly_gain_loss/data/mm10_tss.csv",header = F)

##FOXA1_peak_annotation_peakscore
rownames(w18_all_FOXA1) <- w18_all_FOXA1$PeakID..cmd.annotatePeaks.pl.FOXA1_peak_spike_in_xout_all.mm10.
w18_all_FOXA1 <- w18_all_FOXA1[w18_all_peakscore$X.name..cmd...mergePeaks.over__FOXA1_SKO_spike_in_xout_all.over__FOXA1_DKO_spike_in_xout_all.over__FOXA1_DKO_spike_in_xout_all_FOXA1_SKO_spike_in_xout_all.,]
w18_all_FOXA1$Stat <- w18_all_peakscore$Stat
w18_all_FOXA1 <- w18_all_FOXA1[order(-w18_all_FOXA1$Stat),]

##sc_RNA-seq_DEG_over_FOXA1_ChIP
w18_DEG_0.15_1 <- w18_DEG_0.15[abs(w18_DEG_0.15$avg_log2FC.PF.P.) >= 1,]
w18_DEG_0.15_0.1 <- w18_DEG_0.15[abs(w18_DEG_0.15$avg_log2FC.PF.P.) <= 0.1,]

FOXA1_1.0_DEG_gene <- unique(w18_DEG_0.15_1_FOXA1$Gene[w18_DEG_0.15_1_FOXA1$Gene%in%w18_DEG_0.15_1$Gene])
FOXA1_0.1_unchange_gene <- unique(w18_DEG_0.15_0.1_FOXA1_unchange$Gene[w18_DEG_0.15_0.1_FOXA1_unchange$Gene%in%w18_SKO_only_FOXA1_0.1unchange_peak_highest$Gene.Name])


w18_DEG_0.15_1_FOXA1 <- w18_DEG_0.15_1[w18_DEG_0.15_1$Gene%in%FOXA1_1.0_DEG_gene,]
w18_DEG_0.15_1_FOXA1_up <- w18_DEG_0.15_1_FOXA1[w18_DEG_0.15_1_FOXA1$avg_log2FC.PF.P. > 0,]
w18_DEG_0.15_1_FOXA1_down <- w18_DEG_0.15_1_FOXA1[w18_DEG_0.15_1_FOXA1$avg_log2FC.PF.P. < 0,]
w18_DEG_0.15_0.1_FOXA1_unchange <- w18_DEG_0.15_0.1


w18_SKO_only_FOXA1_1.0up_peak <- w18_all_FOXA1[w18_all_FOXA1$Gene.Name%in%w18_DEG_0.15_1_FOXA1_up$Gene,]
w18_SKO_only_FOXA1_1.0up_peak_highest <- w18_SKO_only_FOXA1_1.0up_peak[!duplicated(w18_SKO_only_FOXA1_1.0up_peak$Gene.Name), ]
w18_SKO_only_FOXA1_1.0down_peak <- w18_all_FOXA1[w18_all_FOXA1$Gene.Name%in%w18_DEG_0.15_1_FOXA1_down$Gene,]
w18_SKO_only_FOXA1_1.0down_peak_highest <- w18_SKO_only_FOXA1_1.0down_peak[!duplicated(w18_SKO_only_FOXA1_1.0down_peak$Gene.Name), ]
w18_SKO_only_FOXA1_0.1unchange_peak <- w18_all_FOXA1[w18_all_FOXA1$Gene.Name%in%w18_DEG_0.15_0.1_FOXA1_unchange$Gene,]
w18_SKO_only_FOXA1_0.1unchange_peak_highest <- w18_SKO_only_FOXA1_0.1unchange_peak[!duplicated(w18_SKO_only_FOXA1_0.1unchange_peak$Gene.Name), ]

##DEGs_not_overlap_FOXA1_binding
xFOXA1_1.0_DEG <- w18_DEG_0.15_1[w18_DEG_0.15_1$Gene%in%FOXA1_1.0_DEG_gene == F,]
xFOXA1_1.0_DEG_up <- xFOXA1_1.0_DEG[xFOXA1_1.0_DEG$avg_log2FC.PF.P. > 0,]
xFOXA1_1.0_DEG_down <- xFOXA1_1.0_DEG[xFOXA1_1.0_DEG$avg_log2FC.PF.P. < 0,]
xFOXA1_0.1_unchange <- w18_DEG_0.15_0.1[w18_DEG_0.15_0.1$Gene%in%FOXA1_0.1_unchange_gene == F,]
xFOXA1_1.0_DEG_up_position_TSS <- TSS_site[TSS_site$V7%in%xFOXA1_1.0_DEG_up$Gene,c(1,2,3,6)]
xFOXA1_1.0_DEG_down_position_TSS <- TSS_site[TSS_site$V7%in%xFOXA1_1.0_DEG_down$Gene,c(1,2,3,6)]
xFOXA1_0.1_unchange_position_TSS <- TSS_site[TSS_site$V7%in%xFOXA1_0.1_unchange$Gene,c(1,2,3,6)]


write.table(w18_SKO_only_FOXA1_1.0up_peak_highest[,c(2,3,4,16)],file = "w18_FOXA1_1.0up_peak_highest.bed", row.names = F, sep = "\t")
write.table(w18_SKO_only_FOXA1_1.0down_peak_highest[,c(2,3,4,16)],file = "w18_FOXA1_1.0down_peak_highest.bed", row.names = F, sep = "\t")
write.table(w18_SKO_only_FOXA1_0.1unchange_peak_highest[,c(2,3,4,16)],file = "w18_SKO_only_FOXA1_0.1unchange_peak_highest.bed", row.names = F, sep = "\t")

write.table(xFOXA1_1.0_DEG_up_position_TSS,file = "xFOXA1_1.0_DEG_up_position_TSS.bed", row.names = F, sep = "\t")
write.table(xFOXA1_1.0_DEG_down_position_TSS,file = "xFOXA1_1.0_DEG_down_position_TSS.bed", row.names = F, sep = "\t")
write.table(xFOXA1_0.1_unchange_position_TSS,file = "xFOXA1_0.1_unchange_position_TSS.bed", row.names = F, sep = "\t")



FOXA1_order_up <- read.csv(file = "~/R_project/FOXA1_luly_gain_loss/data/FOXA1_sorted_up.csv",header = T)
FOXA1_order_unchange <- read.csv(file = "~/R_project/FOXA1_luly_gain_loss/data/FOXA1_sorted_unchange.csv",header = T)
FOXA1_order_down <- read.csv(file = "~/R_project/FOXA1_luly_gain_loss/data/FOXA1_sorted_down.csv",header = T)

rownames(TSS_site) <- TSS_site$V6
FOXA1_order_up_TSS <- TSS_site[FOXA1_order_up$name,]
FOXA1_order_unchange_TSS <- TSS_site[FOXA1_order_unchange$name,]
FOXA1_order_down_TSS <- TSS_site[FOXA1_order_down$name,]

FOXA1_order_up_gene_coord <- gtf_df_gene[gtf_df_gene$gene_name%in%FOXA1_order_up$name,]

FOXA1_order_up_TSS_FOXA1bind <- FOXA1_order_up_TSS[FOXA1_order_up$name%in%w18_SKO_only_FOXA1_1.0up_peak_highest$Gene.Name,]
FOXA1_order_up_TSS_FOXA1bind <- FOXA1_order_up_TSS_FOXA1bind[!duplicated(FOXA1_order_up_TSS_FOXA1bind$V6),]

FOXA1_order_down_TSS_FOXA1bind <- FOXA1_order_down_TSS[FOXA1_order_down$name%in%w18_SKO_only_FOXA1_1.0down_peak_highest$Gene.Name,]
FOXA1_order_down_TSS_FOXA1bind <- FOXA1_order_down_TSS_FOXA1bind[!duplicated(FOXA1_order_down_TSS_FOXA1bind$V6),]

FOXA1_order_unchange_TSS_FOXA1bind <- FOXA1_order_unchange_TSS[FOXA1_order_unchange$name%in%w18_SKO_only_FOXA1_0.1unchange_peak_highest$Gene.Name,]
FOXA1_order_unchange_TSS_FOXA1bind <- FOXA1_order_unchange_TSS_FOXA1bind[!duplicated(FOXA1_order_unchange_TSS_FOXA1bind$V6),]


write.table(FOXA1_order_up_TSS_FOXA1bind[,c(2,3,4)],file = "FOXA1_order_up_TSS_FOXA1bind.bed", row.names = F, sep = "\t")
write.table(FOXA1_order_unchange_TSS_FOXA1bind[,c(2,3,4)],file = "FOXA1_order_unchange_TSS_FOXA1bind.bed", row.names = F, sep = "\t")
write.table(FOXA1_order_down_TSS_FOXA1bind[,c(2,3,4)],file = "FOXA1_order_down_TSS_FOXA1bind.bed", row.names = F, sep = "\t")