library(gridExtra)
library(patchwork)
library(ggplot2)
library(UpSetR)
library(ggpubr)
library(ggsignif)
library(tidyverse)

##Figure 2 ERT epigenomics##

#####
#####processing readcounts and deseq data######

##files available in raw data folder, run either deseq for mll3 (CKO) or for dko
deseq_mll3koha <- read.csv("~/Desktop/data_grhl/epigenomics/diffbindresults_mll3ko_evst.txt_full.txt", sep = '\t', header=TRUE)
readcounts_mll3koha<- readcounts_tss <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_ert_mll3ko_ha.tab")
deseq_oi <- deseq_mll3koha
readcounts_oi <- readcounts_mll3koha

deseq_dkoha <- read.csv("~/Desktop/data_grhl/epigenomics/diffbindresults_dkoha_evst.txt_full.txt", sep = '\t', header=TRUE)
readcounts_dkoha <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_ert_dko_ha.tab", sep = '\t', header=TRUE)
deseq_oi <- deseq_dkoha
readcounts_oi <- readcounts_dkoha


#deseq_oi$chr <- paste0("chr", deseq_oi$seqnames)
#deseq_oi <- cbind(deseq_oi[12], deseq_oi[2:11])
deseq_oi$start <- deseq_oi$start

colnames(readcounts_oi) <- c("chr", "start", "end", 
                              "mll3ko_eth_k4m1", "mll3ko_tam_k4m1","dko_eth_k4m1", "dko_tam_k4m1", 
                              "mll3ko_eth_k4m2", "mll3ko_tam_k4m2","dko_eth_k4m2", "dko_tam_k4m2", 
                              "mll3ko_eth_k4m3", "mll3ko_tam_k4m3","dko_eth_k4m3", "dko_tam_k4m3",
                              "mll3ko_eth_k27a", "mll3ko_tam_k27a","dko_eth_k27a", "dko_tam_k27a", 
                              "mll3ko_eth_h33", "mll3ko_tam_h33","dko_eth_h33", "dko_tam_h33", 
                              "mll3ko_eth_rad21", "mll3ko_tam_rad21","dko_eth_rad21", "dko_tam_rad21",
                              "mll3ko_eth_ha", "mll3ko_tam_ha","dko_eth_ha", "dko_tam_ha", 
                              "mll3ko_eth_p300", "mll3ko_tam_p300","dko_eth_p300", "dko_tam_p300", 
                              "mll3ko_eth_igg", "mll3ko_tam_igg","dko_eth_igg", "dko_tam_igg")
#main wt readcounts processing
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:39] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start

fold_mll3ko_k4m1 <- readcounts_oi2$mll3ko_tam_k4m1 - readcounts_oi2$mll3ko_eth_k4m1
fold_dko_k4m1 <- readcounts_oi2$dko_tam_k4m1 - readcounts_oi2$dko_eth_k4m1
fold_mll3ko_k4m2 <- readcounts_oi2$mll3ko_tam_k4m2 - readcounts_oi2$mll3ko_eth_k4m2
fold_dko_k4m2 <- readcounts_oi2$dko_tam_k4m2 - readcounts_oi2$dko_eth_k4m2
fold_mll3ko_k4m3 <- readcounts_oi2$mll3ko_tam_k4m3 - readcounts_oi2$mll3ko_eth_k4m3
fold_dko_k4m3 <- readcounts_oi2$dko_tam_k4m3 - readcounts_oi2$dko_eth_k4m3

fold_mll3ko_k27a <- readcounts_oi2$mll3ko_tam_k27a - readcounts_oi2$mll3ko_eth_k27a
fold_dko_k27a <- readcounts_oi2$dko_tam_k427a - readcounts_oi2$dko_eth_k27a
fold_mll3ko_h33 <- readcounts_oi2$mll3ko_tam_h33 - readcounts_oi2$mll3ko_eth_h33
fold_dko_h33 <- readcounts_oi2$dko_tam_h33 - readcounts_oi2$dko_eth_h33
fold_mll3ko_rad21 <- readcounts_oi2$mll3ko_tam_rad21 - readcounts_oi2$mll3ko_eth_rad21
fold_dko_rad21 <- readcounts_oi2$dko_tam_rad21 - readcounts_oi2$dko_eth_rad21

fold_mll3ko_ha <- readcounts_oi2$mll3ko_tam_ha - readcounts_oi2$mll3ko_eth_ha
fold_dko_ha <- readcounts_oi2$dko_tam_ha - readcounts_oi2$dko_eth_ha
fold_mll3ko_p300 <- readcounts_oi2$mll3ko_tam_p300 - readcounts_oi2$mll3ko_eth_p300
fold_dko_p300 <- readcounts_oi2$dko_tam_p300 - readcounts_oi2$dko_eth_p300
fold_mll3ko_igg <- readcounts_oi2$mll3ko_tam_igg - readcounts_oi2$mll3ko_eth_igg
fold_dko_igg <- readcounts_oi2$dko_tam_igg - readcounts_oi2$dko_eth_igg

mll3ko_eth_haigg <- readcounts_oi2$mll3ko_eth_ha - readcounts_oi2$mll3ko_eth_igg
mll3ko_tam_haigg <- readcounts_oi2$mll3ko_tam_ha - readcounts_oi2$mll3ko_tam_igg
fold_mll3ko_haigg  <- mll3ko_tam_haigg - mll3ko_eth_haigg

dko_eth_haigg <- readcounts_oi2$dko_eth_ha - readcounts_oi2$dko_eth_igg
dko_tam_haigg <- readcounts_oi2$dko_tam_ha - readcounts_oi2$dko_tam_igg
fold_dko_haigg  <- dko_tam_haigg - dko_eth_haigg

readcounts_oi3 <- cbind(readcounts_oi2, fold_mll3ko_k4m1,fold_dko_k4m1,fold_mll3ko_k4m2,fold_dko_k4m2,
                        fold_mll3ko_k4m3, fold_dko_k4m3, fold_mll3ko_k27a, fold_dko_k27a, fold_mll3ko_h33, fold_dko_h33,
                        fold_mll3ko_rad21, fold_dko_rad21, fold_mll3ko_ha, fold_dko_ha,
                        fold_mll3ko_p300, fold_dko_p300, fold_mll3ko_igg, fold_dko_igg, fold_mll3ko_haigg, fold_dko_haigg)

supermatrix <- left_join(deseq_oi, readcounts_oi3, by = c("seqnames" = "chr", "start", "end"))

#####
#####Overlapping CUT&Tag/CUT&RUN sites with published GRHL2 sites#####

mll3kofull_eo_GRHL2 <- read.csv("~/Desktop/data_grhl/epigenomics/mll3kofull_eo_GRHL2.bed", sep = '\t', header=FALSE)
mll3kofull_eo_GRHL2 <- mll3kofull_eo_GRHL2[1:3]
colnames(mll3kofull_eo_GRHL2) <- c("chr", "start", "end")
supermatrix_grhl2 <- inner_join(supermatrix, mll3kofull_eo_GRHL2, by = c("seqnames" = "chr", "start", "end"))

#subset analysis using grhl2 motitfs from homer output
mll3kofull_eo_GRHL2 <- read.csv("~/Desktop/data_grhl/epigenomics/annotategrhl2_mll3ko_evst.txt", sep = '\t', header=TRUE)
dfmll3kofull_eo_GRHL2 <- mll3kofull_eo_GRHL2[2:4]
colnames(dfmll3kofull_eo_GRHL2) <- c("chr", "start", "end")
dfmll3kofull_eo_GRHL2$start <- dfmll3kofull_eo_GRHL2$start - 1
dfmll3kofull_eo_GRHL2$grhl2motif <- mll3kofull_eo_GRHL2$GRHL2.CP2..HBE.GRHL2.ChIP.Seq.GSE46194..Homer.Distance.From.Peak.sequence.strand.conservation.
dfmll3kofull_eo_GRHL22 <- dfmll3kofull_eo_GRHL2 %>% filter(grhl2motif != "")
supermatrix_grhl2 <- inner_join(supermatrix, dfmll3kofull_eo_GRHL22, by = c("seqnames" = "chr", "start", "end"))

#####
##### significant sites volcano plot #####

supermatrix_grhl2$log10fdr <- -1 * log10(supermatrix_grhl2$FDR)

foldchangedf <- dplyr::select(supermatrix, fold_mll3ko_ha)
foldchangedf$FDR <- supermatrix$FDR
foldchangedf$log10fdr <- -1 * log10(supermatrix$FDR) 
foldchangedf2 <- foldchangedf %>% drop_na()
foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (fold_mll3ko_ha > 1 | fold_mll3ko_ha < -1))

count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_mll3ko_ha > 1 | fold_mll3ko_ha < -1))) #27458 sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_mll3ko_ha > 1))) #13769
count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_mll3ko_ha < -1))) #13689

ggplot() +
  geom_point(data=foldchangedf2, aes(x=fold_mll3ko_ha, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=fold_mll3ko_ha, y=log10fdr), color = "black", size = 0.3) +
  geom_point(data=supermatrix_grhl2, aes(x=fold_mll3ko_ha, y=log10fdr), color = "blue", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,20) +
  xlim(-15,15)


foldchangedf <- dplyr::select(supermatrix, fold_dko_ha)
foldchangedf$FDR <- supermatrix$FDR
foldchangedf$log10fdr <- -1 * log10(supermatrix$FDR) 
foldchangedf2 <- foldchangedf %>% drop_na()
foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (fold_dko_ha > 1 | fold_dko_ha < -1))

count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_dko_ha > 1 | fold_dko_ha < -1))) #27458 sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_dko_ha > 1))) #13769
count(foldchangedf2 %>% filter(FDR < 0.05 & (fold_dko_ha < -1))) #13689

ggplot() +
  geom_point(data=foldchangedf2, aes(x=fold_dko_ha, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=fold_dko_ha, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,20) +
  xlim(-15,15)

#####
##### HA correlation at GRHL2 and ectopic sites #####

readcounts_ertatgrhl2siteschip <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_grhl2sites_cor-3.tab")
readcounts_oi <- readcounts_ertatgrhl2siteschip 

readcounts_ertatgrhl2sites <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_ertatgrhl2sites.tab")
readcounts_oi <- readcounts_ertatgrhl2sites 

readcounts_ectopic <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_ectopic1439sites.tab", sep = '\t', header=FALSE)
readcounts_oi <- readcounts_ectopic

colnames(readcounts_oi) <- c("chr", "start", "end", 
                             "mll3ko_eth_k4m1", "mll3ko_tam_k4m1","dko_eth_k4m1", "dko_tam_k4m1", 
                             "mll3ko_eth_k4m2", "mll3ko_tam_k4m2","dko_eth_k4m2", "dko_tam_k4m2", 
                             "mll3ko_eth_k4m3", "mll3ko_tam_k4m3","dko_eth_k4m3", "dko_tam_k4m3",
                             "mll3ko_eth_k27a", "mll3ko_tam_k27a","dko_eth_k27a", "dko_tam_k27a", 
                             "mll3ko_eth_h33", "mll3ko_tam_h33","dko_eth_h33", "dko_tam_h33", 
                             "mll3ko_eth_rad21", "mll3ko_tam_rad21","dko_eth_rad21", "dko_tam_rad21",
                             "mll3ko_eth_ha", "mll3ko_tam_ha","dko_eth_ha", "dko_tam_ha", 
                             "mll3ko_eth_p300", "mll3ko_tam_p300","dko_eth_p300", "dko_tam_p300", 
                             "mll3ko_eth_igg", "mll3ko_tam_igg","dko_eth_igg", "dko_tam_igg",
                             "chen_grhl2_grhl2ko", "chen_grhl2_wt", "chen_ha_cola1", "chen_ha_grhl2oe",
                             "form_dko_grhl2", "form_wt_grhl2")
#main wt readcounts processing
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:45] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start + 1

fold_mll3ko_k4m1 <- readcounts_oi2$mll3ko_tam_k4m1 - readcounts_oi2$mll3ko_eth_k4m1
fold_dko_k4m1 <- readcounts_oi2$dko_tam_k4m1 - readcounts_oi2$dko_eth_k4m1
fold_mll3ko_k4m2 <- readcounts_oi2$mll3ko_tam_k4m2 - readcounts_oi2$mll3ko_eth_k4m2
fold_dko_k4m2 <- readcounts_oi2$dko_tam_k4m2 - readcounts_oi2$dko_eth_k4m2
fold_mll3ko_k4m3 <- readcounts_oi2$mll3ko_tam_k4m3 - readcounts_oi2$mll3ko_eth_k4m3
fold_dko_k4m3 <- readcounts_oi2$dko_tam_k4m3 - readcounts_oi2$dko_eth_k4m3

fold_mll3ko_k27a <- readcounts_oi2$mll3ko_tam_k27a - readcounts_oi2$mll3ko_eth_k27a
fold_dko_k27a <- readcounts_oi2$dko_tam_k27a - readcounts_oi2$dko_eth_k27a
fold_mll3ko_h33 <- readcounts_oi2$mll3ko_tam_h33 - readcounts_oi2$mll3ko_eth_h33
fold_dko_h33 <- readcounts_oi2$dko_tam_h33 - readcounts_oi2$dko_eth_h33
fold_mll3ko_rad21 <- readcounts_oi2$mll3ko_tam_rad21 - readcounts_oi2$mll3ko_eth_rad21
fold_dko_rad21 <- readcounts_oi2$dko_tam_rad21 - readcounts_oi2$dko_eth_rad21

fold_mll3ko_ha <- readcounts_oi2$mll3ko_tam_ha - readcounts_oi2$mll3ko_eth_ha
fold_dko_ha <- readcounts_oi2$dko_tam_ha - readcounts_oi2$dko_eth_ha
fold_mll3ko_p300 <- readcounts_oi2$mll3ko_tam_p300 - readcounts_oi2$mll3ko_eth_p300
fold_dko_p300 <- readcounts_oi2$dko_tam_p300 - readcounts_oi2$dko_eth_p300
fold_mll3ko_igg <- readcounts_oi2$mll3ko_tam_igg - readcounts_oi2$mll3ko_eth_igg
fold_dko_igg <- readcounts_oi2$dko_tam_igg - readcounts_oi2$dko_eth_igg

mll3ko_eth_haigg <- readcounts_oi2$mll3ko_eth_ha - readcounts_oi2$mll3ko_eth_igg
mll3ko_tam_haigg <- readcounts_oi2$mll3ko_tam_ha - readcounts_oi2$mll3ko_tam_igg
fold_mll3ko_haigg  <- mll3ko_tam_haigg - mll3ko_eth_haigg

dko_eth_haigg <- readcounts_oi2$dko_eth_ha - readcounts_oi2$dko_eth_igg
dko_tam_haigg <- readcounts_oi2$dko_tam_ha - readcounts_oi2$dko_tam_igg
fold_dko_haigg  <- dko_tam_haigg - dko_eth_haigg

readcounts_oi3 <- cbind(readcounts_oi2, fold_mll3ko_k4m1,fold_dko_k4m1,fold_mll3ko_k4m2,fold_dko_k4m2,
                        fold_mll3ko_k4m3, fold_dko_k4m3, fold_mll3ko_k27a, fold_dko_k27a, fold_mll3ko_h33, fold_dko_h33,
                        fold_mll3ko_rad21, fold_dko_rad21, fold_mll3ko_ha, fold_dko_ha,
                        fold_mll3ko_p300, fold_dko_p300, fold_mll3ko_igg, fold_dko_igg, fold_mll3ko_haigg, fold_dko_haigg)
readcounts_oi3ectopic <- cbind(readcounts_oi2, fold_mll3ko_k4m1,fold_dko_k4m1,fold_mll3ko_k4m2,fold_dko_k4m2,
                        fold_mll3ko_k4m3, fold_dko_k4m3, fold_mll3ko_k27a, fold_dko_k27a, fold_mll3ko_h33, fold_dko_h33,
                        fold_mll3ko_rad21, fold_dko_rad21, fold_mll3ko_ha, fold_dko_ha,
                        fold_mll3ko_p300, fold_dko_p300, fold_mll3ko_igg, fold_dko_igg, fold_mll3ko_haigg, fold_dko_haigg)

ggplot() +
  geom_abline(slope = 1) +
  #geom_point(data=readcounts_oi3ectopic, aes(x=mll3ko_tam_ha, y=dko_tam_ha), color = "gray55", size = 0.4) +
  geom_point(data=readcounts_oi3, aes(x=mll3ko_tam_ha, y=dko_tam_ha), color = "blue", size = 0.4) +
  #geom_point(data=readcounts_oi3, aes(x=dko_tam_igg, y=dko_tam_ha), color = "grey", size = 0.4) +
  ylim(c(0,20)) +
  xlim(c(0,20)) +
theme_bw() +
  theme(aspect.ratio = 1)

ct <- cor.test(readcounts_oi3$mll3ko_tam_ha, readcounts_oi3$dko_tam_ha, method = "pearson") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)
ct <- cor.test(readcounts_oi3ectopic$mll3ko_tam_ha, readcounts_oi3ectopic$dko_tam_ha, method = "pearson") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)

ggplot() +
  geom_abline(slope = 1) +
  #geom_point(data=readcounts_oi3ectopic, aes(x=mll3ko_eth_ha, y=dko_eth_ha), color = "gray55", size = 0.4) +
  geom_point(data=readcounts_oi3, aes(x=mll3ko_eth_ha, y=dko_eth_ha), color = "blue", size = 0.4) +
  ylim(c(0,20)) +
  xlim(c(0,20)) +
  theme_bw() +
  theme(aspect.ratio = 1)

count(readcounts_oi3 %>% filter(mll3ko_eth_ha > 12 & dko_eth_ha > 12))
count(readcounts_oi3 %>% filter(mll3ko_tam_ha > 12 & dko_tam_ha > 12))

#readcounts_oi3v2 <- readcounts_oi3 %>% filter(mll3ko_tam_ha > 1)

ggplot() +
  #geom_abline(slope = 1) +
  #geom_point(data=readcounts_oi3ectopic, aes(x=mll3ko_eth_ha, y=dko_eth_ha), color = "gray55", size = 0.4) +
  geom_point(data=readcounts_oi3v2, aes(x=chen_grhl2_wt, y=chen_ha_grhl2oe), color = "blue", size = 0.4) +
  ylim(c(7, 15)) +
  xlim(c(7,15)) +
  theme_bw() +
  theme(aspect.ratio = 1)

ct <- cor.test(readcounts_oi3v2$mll3ko_tam_ha, readcounts_oi3v2$chen_ha_grhl2oe, method = "spearman") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)

ct <- cor.test(readcounts_oi3v2$chen_grhl2_wt, readcounts_oi3v2$chen_ha_grhl2oe, method = "pearson") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)

ct <- cor.test(readcounts_oi3$chen_grhl2_wt, readcounts_oi3$form_wt_grhl2, method = "spearman") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)

#####
#####creating venn diagrams for filtering#####

library(VennDiagram)

grid.newpage()                                        
draw.pairwise.venn(5135,3013,1483)

grid.newpage()                                        
draw.pairwise.venn(5135,332,110)

grid.newpage()                                        
draw.pairwise.venn(3013,332,146)


#####
##### Filtering and numbers for background enhancer set #####

mll3kosig_k27a_homer <- read.delim("~/Desktop/data_grhl/epigenomics/k27a_full.homer.txt")
homerpeakannotationcleanup <- function(mergedannottable) {
  
  # df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "Other", x)))
  
  df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "TTS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "5'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "5'UTR", x)))
  
  #df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("p.*", "Promoter-TSS", x)))
  # 
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_1$", "cluster_01", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_2", "cluster_02", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_3", "cluster_03", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_4", "cluster_04", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_5", "cluster_05", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_6", "cluster_06", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_7", "cluster_07", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_8", "cluster_08", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_9", "cluster_09", x)))
}
mll3kosig_k27a_homer2 <- homerpeakannotationcleanup(mll3kosig_k27a_homer)
mll3kosig_k27a_homer2.ii <- mll3kosig_k27a_homer2 %>% filter((Annotation == "Intron" | Annotation == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
#write.table(mll3kosig_k27a_homer2.ii[2:4], "~/Desktop/mll3kosig_k27a_homer2.ii.bed", sep='\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

#####
#####genome feature annotation#####
homerpeakannotationcleanup <- function(mergedannottable) {
  
  df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "TTS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "5'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "5'UTR", x)))
}
x <- read.delim("~/Desktop/data_grhl/epigenomics/grhl2_332peaks.sort.homer.bed")

homerpeakannotationcleanup <- function(mergedannottable) {
  
  # df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "Other", x)))
  
  df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "TTS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "5'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "5'UTR", x)))
  
  #df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("p.*", "Promoter-TSS", x)))
  # 
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_1$", "cluster_01", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_2", "cluster_02", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_3", "cluster_03", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_4", "cluster_04", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_5", "cluster_05", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_6", "cluster_06", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_7", "cluster_07", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_8", "cluster_08", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_9", "cluster_09", x)))
}
x <- homerpeakannotationcleanup(x)
x <- as.data.frame(x[, "Annotation"])
colnames(x) <- "Annotation"
x$peakset <- "endogenous"

xx <- read.delim("~/Desktop/data_grhl/epigenomics/eo_sigup_merge.notgrhl2.homer.txt")
xx <- homerpeakannotationcleanup(xx)
xx <- as.data.frame(xx[, "Annotation"])
colnames(xx) <- "Annotation"
xx$peakset <- "ectopic"

xxx <- rbind(x, xx)

ggplot(xxx, aes(x=peakset, fill = Annotation)) +
  geom_bar(width = 0.9, position="fill") +
  theme_classic() +
  theme(text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))



#####








