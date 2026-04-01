library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(patchwork)
library(VennDiagram)


#####Import supermatrix data from Boileau 2023#####

supermatrix <- read.delim("~/Desktop/data/processed_data/supermatrix_k27a.20220427.txt", sep = '\t', header = TRUE)

## use with wt data
k27a_atac_filter <- read.table("~/Desktop/k27a_at_peaksfiltered.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

#distal wt only
supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

##mll34 dependent or not
supermatrix.up.k27adown <- supermatrix.up %>% filter(pk_fold_f_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.up.k27anon <- supermatrix.up %>% filter(pk_fold_f_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.non.k27adown <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.non.k27anon <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

#####
#####Generate deeptools compatible bed file#####

generate_deepheatmapgrouped <- function(listofinputs, outputfilename){
  "input a list of bed file coordinates chr,start,end and output a single txt file with # separating groups for use with deeptools plot heatmap"
  stuffer <- "#"
  write.table(listofinputs[[1]][1:3], outputfilename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  for(i in 2:length(listofinputs)){
    nextgroup <- listofinputs[[i]][1:3]
    write.table(stuffer, outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
    write.table(nextgroup, outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t') 
    
  }
}

generate_deepheatmapgrouped(list(supermatrix.non.k27anon, supermatrix.non.k27adown, supermatrix.up.k27anon, supermatrix.up.k27adown), "~/Desktop/data_grhl/epigenomics/ms2_formativek27a.bed")

#####
#####Generate gimme compatible df#####

sortedbed <- read.delim("~/Desktop/data_grhl/epigenomics/motifcounts/RBRB07CT1_k27a_cpm.heatmap.up.1130550.bed")
sortedbed_g1 <- sortedbed %>% filter(deepTools_group == "yy4_use.bed")
sortedbed_g2 <- sortedbed %>% filter(deepTools_group == "genes")

y_one <- sortedbed_g1
colnames(y_one) <- c("chr", "start", "end")
y_one$peaktype <- "up.k27anon"
y_two <- sortedbed_g2
colnames(y_two) <- c("chr", "start", "end")
y_two$peaktype <- "up.k27adown"

x <- rbind(y_one, y_two) 
gimme <- as.data.frame(paste0(x$chr, ":", x$start, "-", x$end))
gimme$cluster <- x$peaktype
colnames(gimme) <- c("loc", "cluster")
#write.table(gimme, "~/Desktop/gimme2024y4peaksv2.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')

#####
#####Gimmemotifs maelstrom output#####

#gimmeoutput <- read.delim("~/Desktop/data_grhl/epigenomics/final.out.4hocomoco.txt")
gimmeoutput <- read.delim("~/Desktop/data_grhl/epigenomics/final.out.4homer.txt")
gimmeoutput <- read.delim("~/Desktop/data_grhl/epigenomics/final.out.4jaspar.txt")
gimmeoutput <- read.delim("~/Desktop/data_grhl/epigenomics/final.out.4image.txt")

colnames(gimmeoutput) <- c("motif", "zscore.y4ind", "zscore.y4dep", 
                           "percent.y4dep", "percent.y4ind")
rownames(gimmeoutput) <- gimmeoutput$motif
gimmeoutput2 <- gimmeoutput %>% filter(percent.y4ind > 5 | percent.y4dep > 5)

col_fun1 = colorRamp2(c(0, 15), c("white", "magenta"))
col_fun2 = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

ht1 <- Heatmap(gimmeoutput2[4:5], col_fun1, 
               #cluster_columns = FALSE,
               row_order = order(-gimmeoutput2$percent.y4dep),
               show_row_names = TRUE,
               row_names_side = "left",
               column_order = sort(c("percent.y4ind", "percent.y4dep")),
               row_names_gp = gpar(fontsize = 8),
               height = nrow(gimmeoutput2[4:5])*unit(5, "mm"),
               width = ncol(gimmeoutput2[4:5])*unit(10, "mm"))

ht2 <- Heatmap(gimmeoutput2[2:3], col_fun2,
               show_row_names = FALSE,
               #cluster_columns = FALSE,
               column_order = sort(c("zscore.y4ind", "zscore.y4dep")),
               height = nrow(gimmeoutput2[4:5])*unit(5, "mm"),
               width = ncol(gimmeoutput2[4:5])*unit(10, "mm"))
ht1 + ht2






#####




#####tie gimmeoutput only y4grid######

motifcount <- read.delim("~/Desktop/data_grhl/epigenomics/motifcounts/motif.count.4homer.txt")
motifcount2 <- as.data.frame(c(motifcount["loc"], 
                               motifcount["ATGCATWATGCATRW_OCT.OCT.short.POU.Homeobox..NPC.OCT6.ChIP.Seq.GSE43916..Homer"],
                               motifcount["AAACYKGTTWDACMRGTTTB_GRHL2.CP2..HBE.GRHL2.ChIP.Seq.GSE46194..Homer"],
                               motifcount["ACWTCAAAGG_TCFL2.HMG..K562.TCF7L2.ChIP.Seq.GSE29196..Homer"],
                               motifcount["AAGGKGRCGCAGGCA_ZNF165.Zf..WHIM12.ZNF165.ChIP.Seq.GSE65937..Homer"],
                               motifcount["ATTTGCATAACAATG_OCT4.SOX2.TCF.NANOG.POU.Homeobox.HMG..mES.Oct4.ChIP.Seq.GSE11431..Homer"],
                               motifcount["NDCTAATTAS_En1.Homeobox..SUM149.EN1.ChIP.Seq.GSE120957..Homer"],
                               motifcount["TGCTGAGTCA_Bach2.bZIP..OCILy7.Bach2.ChIP.Seq.GSE44420..Homer"]))
colnames(motifcount2) <- c("loc", "Oct6", "GRHL2", "TCF7L2", "ZNF165", "OSTN", "EN1", "Bach2")
motifcountx <- data.frame(do.call('rbind', strsplit(as.character(motifcount2$loc),':',fixed=TRUE)))
motifcountxx <- data.frame(do.call('rbind', strsplit(as.character(motifcountx$X2),'-',fixed=TRUE)))
motifcountxxx <- as.data.frame(c(motifcountx["X1"], motifcountxx["X1"], motifcountxx["X2"]))
colnames(motifcountxxx) <- c("chr", "start", "end")
motifcountxxx$start <- as.numeric(motifcountxxx$start)
motifcountxxx$end <- as.numeric(motifcountxxx$end)
motifcount3 <- cbind(motifcountxxx, motifcount2)
motifmaster <- motifcount3


sortedbed$rank <- seq(1,nrow(sortedbed))
sortedbed_combined <- inner_join(sortedbed, motifmaster, by = c("X.chrom" = "chr", "start" = "start", "end" = "end"))

sortedbed_combined$grouping <- c(rep(0,271), rep(1,901))


col_fun1 = colorRamp2(c(0, 1, 2), c("gray", "cyan", "blue"))
col_fun1 = colorRamp2(c(1, 2, 3), c("white", "cyan", "blue"))
col_fun1 = colorRamp2(c(0, 1, 3), c("gray85", "blue", "black"))
col_fun1 = colorRamp2(c(0, 1, 3), c("gray95", "black", "black"))

ht1 <- Heatmap(sortedbed_combined[16:23], col_fun1, 
               cluster_columns = FALSE,
               row_order = order(sortedbed_combined$rank),
               show_row_names = FALSE,
               border_gp = gpar(col = "black", lty = 1),
               width = ncol(sortedbed_combined[16:23])*unit(10, "mm"))
ht1







#####
#####rna-seq lcpm scatterplots#####

#starting from filteredlcpm boileau2023 rnaseq
genesofinterest <- c("Grhl2")
genesofinterest <- c("Grhl2","Pou3f1", "Tcf7l2", "Znf165", "Pou5f1", "En1", "Bach2")

filteredlcpmhighlight <- subset(filteredlcpm, rownames(filteredlcpm) %in% genesofinterest)
filteredlcpmhighlight$genesymbol <- rownames(filteredlcpmhighlight)

###form lcpm plot
compare_group_first=c("FWT", "FWT", "FWT", "FCKO")
compare_group_firstname=c("Form WT", "Form WT","Form WT", "Form MLL3KO")
compare_group_second=c("FCKO", "FDKO", "FdCD", "FDKO")
compare_group_secondname=c("Form MLL3KO", "Form DKO", "Form dCD", "Form MLL3KO")
compare_group_z <- c("FCKO_WT", "FDKO_WT", "FdCD_WT", "FDKO_FCKO")

siggene_group <- list(k4, k5, k6, k10)
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname, compare_group_z))
siggene_group_counter = 0
plotlist=list()
combinedplotlist=list()

plotarray <- for(i in 1:nrow(comparematrix)) {
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  zgroup <- as.character(row[1,5])
  siggene_group_counter = siggene_group_counter + 1

  siggenes <- merge(siggene_group[siggene_group_counter],filteredlcpm, by = 'ensembllistanno')
  #zgenes <- eval(parse(text = zgroup))
  filteredlcpm_wtsig.up <- siggenes %>% filter(eval(parse(text = zgroup)) > 1)
  filteredlcpm_wtsig.down <- siggenes %>% filter(eval(parse(text = zgroup)) < -1)
  
  countup <- filteredlcpm_wtsig.up %>% count(eval(parse(text = zgroup)) > 1)
  countup <- paste(countup[1,2], "Up", sep=' ')
  grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
                              gp=gpar(col="magenta3", fontsize=13)))
  
  #countdown <- filteredlcpm_wtsig.down %>% count(rnafold_geno_wt < -1)
  
  countdown <- filteredlcpm_wtsig.down %>% count(eval(parse(text = zgroup)) < -1)
  countdown <- paste(countdown[1,2], "Down", sep=' ')
  grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
                                gp=gpar(col="magenta3", fontsize=13)))
  
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot(data=filteredlcpm, aes_string(x=xsample, y=ysample)) +
                                                  geom_point(size=0.2, shape=23, color = "gray") +
                                                  geom_point(data=filteredlcpm_wtsig.up, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.down, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.2) +
                                                  geom_text_repel(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample, label="genesymbol", size = 12), color = "black", force = 175, box.padding = 1, max.overlaps = 100)+ #, color = "black", force = 30) +
                                                  geom_point(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample), color = 'black') +
                                                  annotation_custom(grobup) +
                                                  annotation_custom(grobdown) +
                                                  xlab(paste(xname, "(Log2 CPM)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  combinedplotlist[[siggene_group_counter]] <- plotlist[[siggene_group_counter]]
  
}
wrap_plots(combinedplotlist)





####

#####
#####grhl2 sites readcounts scatterplots#####

readcounts_nf_grhl2sites <- read.table('~/Desktop/data_grhl/epigenomics/readcounts/readcounts_grhl2sites.tab', sep = '\t', header=TRUE)
readcounts_oi <- readcounts_nf_grhl2sites
  
#main grhl2 site readcounts labels
colnames(readcounts_oi) <- c("chr", "start", "end", "N_WT_Grhl2", "F_WT_Grhl2", "N_DKO_Grhl2", "F_DKO_Grhl2",
                             "N_WT_k4m1", "F_WT_k4m1","N_CKO_k4m1", "F_CKO_k4m1", "N_DKO_k4m1", "F_DKO_k4m1",
                             "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a",
                               "CR2_IgG", "N_WT_igg", "F_WT_igg", "F_WT_MLL4")

#convert MLL4 dko sub negative values into 0s
readcounts_oi$F_WT_MLL4 <- pmax(readcounts_oi$F_WT_MLL4, 0)
#main wt readcounts processing
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:21] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start + 1

ggplot(data=readcounts_oi2, aes(x=F_WT_Grhl2, y=F_DKO_Grhl2)) +
  geom_point(size = 0.2) +
  geom_smooth(method="lm") +
  xlim(c(0,15)) +
  ylim(c(0,15)) +
  theme_bw()+
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))
  
ct <- cor.test(readcounts_oi2$F_WT_Grhl2, readcounts_oi2$F_DKO_Grhl2, method = "pearson") #, exact = FALSE
ct_pearson <- round(as.numeric(ct[4]), digits = 2)


#####
#####GRHL2 cutandrun diffbind #####

deseq_oi <- read.csv("~/Desktop/data_grhl/epigenomics/diffbindresults_nf_WT_GRHL2.txt_full.txt", sep = '\t', header=TRUE)
readcounts_oi <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_grhl2counts.tab", sep = '\t', header=FALSE)
colnames(readcounts_oi) <- c("chr", "start", "end", 
                             "Naive_WT_Grhl2", "Form_WT_Grhl2", "Naive_DKO_Grhl2", "Form_DKO_Grhl2")
#main wt readcounts processing
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:7] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start

fold_wtgrhl2 <- readcounts_oi2$Form_WT_Grhl2 - readcounts_oi2$Naive_WT_Grhl2
readcounts_oi3 <- cbind(readcounts_oi2, fold_wtgrhl2)
supermatrix <- left_join(deseq_oi, readcounts_oi3, by = c("seqnames" = "chr", "start", "end"))

significant <- supermatrix %>% filter(FDR < 0.05 & abs(fold_wtgrhl2) > 1)
count(supermatrix %>% filter(FDR < 0.05 & fold_wtgrhl2 > 1)) #391
count(supermatrix %>% filter(FDR < 0.05 & fold_wtgrhl2 < -1)) #2949
nrow(supermatrix) #131329
ggplot() +
  geom_point(data=supermatrix, aes(x=fold_wtgrhl2, y=-log10(FDR)), color = "gray", size = 0.3)+
  geom_point(data=significant, aes(x=fold_wtgrhl2, y=-log10(FDR)), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,15) +
  xlim(-10,10)





#####

