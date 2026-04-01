library(limma)
library(DESeq2)
library(edgeR)
library(stringr)
library(dplyr)
library(ggplot2)
library(GGally)
library(patchwork)
library(gridExtra)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(UpSetR)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(umap)
library(VennDiagram)
#raw input files avaible in raw data folder
diffgenes_cko <- read.table("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_16hr_ckotam_vs_ckoeth.20230806.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                            stringsAsFactor=FALSE)
diffgenes_cko$gene.name <- rownames(diffgenes_cko)

#####Generate log2cpm TMM normalized count sheet######
# read data
#use
counts8hr <- read.table("~/Desktop/data_grhl/grhl2_ertrnaseq/salmon.merged.gene_counts_r7v2.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                        stringsAsFactor=FALSE)
colnames(counts8hr) <- paste0("8hr_",colnames(counts8hr))
counts16hr <- read.table("~/Desktop/data_grhl/grhl2_ertrnaseq/salmon.merged.gene_counts_r8v2.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                         stringsAsFactor=FALSE)
colnames(counts16hr) <- paste0("16hr_",colnames(counts16hr))
xx <- full_join(counts8hr, counts16hr, by = c("8hr_gene_name" = "16hr_gene_name"))
counts <- xx[2:31]
rownames(counts) <- gene_name <- xx$`8hr_gene_name`
#output <- "~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenesdkoeth_vs_dkotam.20230427.anno.csv"

#colnames(counts) <- c("A1", "A2", "A3", "B1","B2", "B3")
construct <- c("A", "A", "A", "B", "B", "B", "C", "C", "C","D", "D", "D", "E", "E","E")
construct <- c(paste0(construct,"16hr"), paste0(construct,"8hr"))
replicate <- c("1", "2", "3")
colnames(counts) <- c(paste0(construct,replicate))
annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
rownames(annot) <- c(paste0(construct,replicate))
####
colnames(counts) <- c(rownames(annot))

y <- DGEList(counts=counts,samples=annot)
samplenames <- colnames(y) 
samplenames

construct <- as.factor(y$samples$construct)
#time <- as.factor(y$samples$time)
replicate <- as.factor(y$samples$replicate)
y$geneid <- rownames(y)
lib.size <- as.factor(y$samples$lib.size)

# Subset samples:
x <- y[,,keep.lib.sizes=TRUE] ##subset by what?
dim(x)

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

table(rowSums(x$counts==0)==6)

keep.exprs <- rowSums(cpm>1)>=6 ### changed from at least 1 cpm in 2 samples to
x <- x[keep.exprs,, keep.lib.sizes=FALSE] # (=F) recounts lib size after trim

dim(x)

## Plot Log-CPM RAW and Filtered data
nsamples <- ncol(x)
col <- brewer.pal(8, "Paired") # find pallete to accommodate >12 samples if necessary
par(mfrow=c(1,2))

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

lcpm_x <- cpm(x, log=TRUE)
x_norm <- calcNormFactors(x, method = "TMM")
x_norm$samples$norm.factors

lcpm_x_norm <- cpm(x_norm, log=TRUE)
x_norm <- x
x_norm$samples$norm.factors <- 1

par(mfrow=c(1,2))
lcpm <- cpm(x_norm, log=TRUE)

## Boxplots of Data pre and post-normalization
boxplot(lcpm_x, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
x_norm <- calcNormFactors(x_norm)
x_norm$samples$norm.factors

lcpm <- cpm(x, log=TRUE)
boxplot(lcpm_x_norm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

df_lcpm <- as.data.frame(lcpm)

#batch correction via limma
#batch <- c(rep("A", 15),rep("B", 15))
#df_lcpmbatch <- as.data.frame(removeBatchEffect(df_lcpm, batch=batch))
#df_lcpm2 <- df_lcpmbatch
#skip batch correction
df_lcpm2 <- df_lcpm

#process fold changes and add to dataframe
df_lcpm2 <- mutate(df_lcpm2, mll3koeth = rowMeans(dplyr::select(df_lcpm2, starts_with("A8")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, mll3kotam = rowMeans(dplyr::select(df_lcpm2, starts_with("B8")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, mll3kowt = rowMeans(dplyr::select(df_lcpm2, starts_with("C8")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, dkoeth = rowMeans(dplyr::select(df_lcpm2, starts_with("D8")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, dkotam = rowMeans(dplyr::select(df_lcpm2, starts_with("E8")), na.rm = TRUE))

df_lcpm2 <- mutate(df_lcpm2, mll3koeth16 = rowMeans(dplyr::select(df_lcpm2, starts_with("A16")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, mll3kotam16 = rowMeans(dplyr::select(df_lcpm2, starts_with("B16")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, mll3kowt16 = rowMeans(dplyr::select(df_lcpm2, starts_with("C16")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, dkoeth16 = rowMeans(dplyr::select(df_lcpm2, starts_with("D16")), na.rm = TRUE))
df_lcpm2 <- mutate(df_lcpm2, dkotam16 = rowMeans(dplyr::select(df_lcpm2, starts_with("E16")), na.rm = TRUE))

df_lcpm2$fold_mll3ko_tameth <- df_lcpm2$mll3kotam - df_lcpm2$mll3koeth  
df_lcpm2$fold_dko_tameth <- df_lcpm2$dkotam - df_lcpm2$dkoeth 
df_lcpm2$fold_mll3ko_wteth <- df_lcpm2$mll3kowt - df_lcpm2$mll3koeth 

df_lcpm2$fold_mll3ko_tameth16 <- df_lcpm2$mll3kotam16 - df_lcpm2$mll3koeth16
df_lcpm2$fold_dko_tameth16 <- df_lcpm2$dkotam16 - df_lcpm2$dkoeth16 
df_lcpm2$fold_mll3ko_wteth16 <- df_lcpm2$mll3kowt16 - df_lcpm2$mll3koeth16 

df_lcpm_ert <- df_lcpm2
df_lcpm_ert$gene.name <- rownames(df_lcpm_ert)
#write.csv(df_lcpm_ert, "~/Desktop/data_grhl/grhl2_ertrnaseq/df_lcpm_ert_20240421.csv")

#####
#####Differential gene expression using DESeq2#####

counts <- read.table("~/Desktop/data_grhl/grhl2_ertrnaseq/salmon.merged.gene_counts_r7v2.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- read.table("~/Desktop/data_grhl/grhl2_ertrnaseq/salmon.merged.gene_counts_r8v2.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- counts[2:16]
ckoeth <- counts[, 1:3]
ckotam <- counts[, 4:6]
ckowt <- counts[, 7:9]
dkoeth <- counts[, 10:12]
dkotam <- counts[, 13:15]
#dkotam <- counts[,c("dko_tam_REP1", "dko_tam_REP3")]


#edit counts here to change which pair-wise differential gene expression to perform
counts <- data.frame(dkoeth, dkotam)
counts <- data.frame(ckoeth, ckowt)
output <- "~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_8hr_ckowt_vs_ckoeth.20240407.anno.csv"
output <- "~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_8hr_dkotam_vs_dkoeth.20240421.anno.csv"
colnames(counts) <- c("A1", "A2", "A3", "B1","B2", "B3")
construct <- c("A", "A", "A", "B", "B", "B")

replicate <- c("1", "2", "3")
annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
rownames(annot) <- c("A1","A2","A3","B1","B2","B3")
#rownames(annot) <- c("A1","A2","A3","B1","B2")
DESeqOurcts <- as.matrix(counts)
DESeqOurcts <- apply (DESeqOurcts, c (1, 2), function (x) {
  (as.integer(x)) })

DESeqOurcoldata <- as.matrix(annot)

all(rownames(DESeqOurcoldata) %in% colnames(DESeqOurcts))
all(rownames(DESeqOurcoldata) == colnames(DESeqOurcts))

DESeqOurdds <- DESeqDataSetFromMatrix(countData = DESeqOurcts,
                                      colData = DESeqOurcoldata,
                                      design = ~ construct)
DESeqOurdds
# add additional metadata
DESeqOurfeatureData <- data.frame(gene=rownames(DESeqOurcts))
mcols(DESeqOurdds) <- DataFrame(mcols(DESeqOurdds), DESeqOurfeatureData)
mcols(DESeqOurdds)
# prefiltering
DESeqOurdds <- DESeqOurdds[ rowSums(counts(DESeqOurdds)) > 1, ]
# explicitly set reference (default is 1st by alphabetic); not necessary if using contrasts
##DESeqOurdds$condition <- relevel(DESeqOurdds$group, ref="e75") ##obvs 7.5 is not in this but what would be good ref?
### Differential Expression
OurddsMF <- DESeqOurdds

design(OurddsMF) <- formula(~ construct)
OurddsMF <- DESeq(OurddsMF)
OurresMF <- results(OurddsMF)
head(OurresMF)

# MA plot
plotMA(OurresMF, alpha = 0.05, ylim=c(-8,8))

resMFOrdered <- OurresMF[order(OurresMF$padj),] # reordered by padj
summary(resMFOrdered)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05

#filter reads for padj and log2FC and generate .csv file named by output above

resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))
resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05 & (log2FoldChange > 1))
resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05)

write.csv(resMFOrderedfiltered, output)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05
#converted <- convertensemblgenetosymbol(resMFOrderedfiltered)
#write.csv(converted, output)

#####
#####skip upstream processing, load ert dfs#####

df_lcpm_ert <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/df_lcpm_ert_20240421.csv", sep = ',', row.names = 1)
diffgeneset_ckowt8hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_8hr_ckowt_vs_ckoeth.20240407.anno.csv")
diffgeneset_cko8hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_8hr_ckotam_vs_ckoeth.20240308.anno.csv")
diffgeneset_dko8hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_8hr_dkotam_vs_dkoeth.20240308.anno.csv")

diffgeneset_ckowt16hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_16hr_ckowt_vs_ckoeth.20240312.anno.csv")
diffgeneset_cko16hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_16hr_ckotam_vs_ckoeth.20230806.anno.csv")
diffgeneset_dko16hr <- read.csv("~/Desktop/data_grhl/grhl2_ertrnaseq/Diffgenes_16hr_dkotam_vs_dkoeth.20240308.anno.csv")

#####
##### scatterplots for 8hr ert #####
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt8hr$X & fold_mll3ko_wteth < -0.5)) #6 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt8hr$X & fold_mll3ko_wteth > 0.5)) #23 up
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth < -0.5)) #15 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0.5)) #101 up
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & fold_dko_tameth < -0.5)) #5 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & fold_dko_tameth > 0.5)) #9 up

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt8hr$X & abs(fold_mll3ko_wteth) > 0.5)
genestohighlight <- c("Kmt2d", "Grhl2")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=mll3koeth, y=mll3kowt), size = 0.2, color = "gray") +
  geom_point(data=genesubset , aes(x=mll3koeth, y=mll3kowt), size = 0.4, color = "orange2") +
  geom_text_repel(data = highlight, aes(x=mll3koeth, y=mll3kowt, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +   
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")



genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & abs(fold_mll3ko_tameth) > 0.5)
genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Wnt7b", "Cldn4", "Epcam")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=mll3koeth, y=mll3kotam), size = 0.2, color = "gray") +
  geom_point(data=genesubset , aes(x=mll3koeth, y=mll3kotam), size = 0.4, color = "orange2") +
  geom_text_repel(data = highlight, aes(x=mll3koeth, y=mll3kotam, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & abs(fold_dko_tameth) > 0.5)
genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Cldn4", "Wnt7b", "Epcam")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

genesubset2 <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & abs(fold_dko_tameth) > 0.5)


ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=dkoeth, y=dkotam), size = 0.2, color = "gray") +
  geom_point(data=genesubset , aes(x=dkoeth, y=dkotam), size = 0.4, color = "orange2") +
  #geom_point(data=genesubset3 , aes(x=dkoeth, y=dkotam), size = 0.4, color = "orange2") +
  geom_text_repel(data = highlight, aes(x=dkoeth, y=dkotam, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

#####
#####UMAP analysis on samples ######

#mll3koert genes only, choose either 8hr or 18hr
xx <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X |
                               gene.name %in% diffgeneset_dko8hr$X)
xx <- xx[16:30]
 
# xx <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X |
#                                gene.name %in% diffgeneset_dko16hr$X)
# xx <- xx[1:15]

filteredreplicatelcpmt <- as.data.frame(t(xx))
counts.umap <- umap(filteredreplicatelcpmt, n_neighbors = 4) #n_neighbors = 5
head(counts.umap$layout)
sample <- c(rep("mll3koeth",3), rep("mll3kotam",3), rep("mll3kowt",3), rep("dkoeth",3), rep("dkotam",3))

replicatename <- sample

layoutdf <- as.data.frame(counts.umap$layout)
colnames(layoutdf) <- c("one", "two")
layoutdf$samplename <- replicatename

umap_8h <- ggplot(layoutdf, mapping = aes(x=one , y=two,fill = samplename)) +
  scale_shape_manual(values = c(21, 21, 22, 22, 23)) +
  geom_point(size=3, aes(shape=samplename)) +
  theme_bw() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2))
umap_8h


#####
##### scatterplots for 16hr ert #####

nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt16hr$X & fold_mll3ko_wteth16 < -0.5)) #5 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt16hr$X & fold_mll3ko_wteth16 > 0.5)) #4 up
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 < -0.5)) #42 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 > 0.5)) #135 up
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & fold_dko_tameth16 < -0.5)) #0 down
nrow(df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & fold_dko_tameth16 > 0.5)) #3 up

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt16hr$X & abs(fold_mll3ko_wteth16) > 0.5)
genestohighlight <- c("Kmt2d", "Grhl2")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight)

ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=mll3koeth16, y=mll3kowt16), size = 0.2, color = "gray") +
  geom_point(data=genesubset , aes(x=mll3koeth16, y=mll3kowt16), size = 0.4, color = "orange2") +
  geom_text_repel(data = highlight, aes(x=mll3koeth16, y=mll3kowt16, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +   
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")



genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & abs(fold_mll3ko_tameth16) > 0.5)
genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Cdh1", "Wnt7b", "Cldn4", "Epcam")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=mll3koeth16, y=mll3kotam16), size = 0.2, color = "gray") +
  geom_point(data=genesubset , aes(x=mll3koeth16, y=mll3kotam16), size = 0.4, color = "orange2") +
  geom_text_repel(data = highlight, aes(x=mll3koeth16, y=mll3kotam16, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & abs(fold_dko_tameth16) > 0.5)
genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Cldn4", "Wnt7b", "Epcam", "Cdh1")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))
genesubset2 <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & abs(fold_dko_tameth16) > 0.5)

ggplot()+
  geom_point(data=df_lcpm_ert, aes(x=dkoeth16, y=dkotam16), size = 0.2, color = "gray") +
  geom_point(data=genesubset2 , aes(x=dkoeth16, y=dkotam16), size = 0.4, color = "orange") +
  geom_text_repel(data = highlight, aes(x=dkoeth16, y=dkotam16, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  ylim(c(-3, 15)) +
  xlim(c(-3, 15)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")


#####
#####boxplots for mll3ko ert 8hr sig #####
genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0.5)
# genesubset_ert_targets <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & gene.name %in% geneset_targets$ensembllistanno)
# genesubset_ert_targetsxformup <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & gene.name %in% geneset_targetsxformup$ensembllistanno)
# genesubset_ert_nearest <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & gene.name %in% geneset_nearest$ensembllistanno)
# genesubset_ert_targetsxnearest <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & gene.name %in% geneset_targetsxnearest$ensembllistanno)
genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 > 0.5)


all <- df_lcpm_ert
all$type <- "all"
subset1 <- genesubset
subset1$type <- "ERT"
# subset2 <- genesubset_ert_targets
# subset2 <- genesubset_ert_nearest
# subset2$type <- "targets"


y0 <- gather(all[,c("fold_mll3ko_tameth", "fold_dko_tameth")])
y0$type <- "all"
y1 <- gather(subset1[,c("fold_mll3ko_tameth", "fold_dko_tameth")])
y1$type <- "subset1"
# y2 <- gather(subset2[,c("fold_mll3ko_tameth", "fold_dko_tameth")])
# y2$type <- "subset2"

z <- rbind(y0, y1)
z$keytype <- paste0(z$type, z$key)

#z$key <- factor(z$key, levels = c("FCKO_WT", "FDKO_WT"))
#z$keytype <- factor(z$keytype, levels = c("allfold_mll3ko_tameth", "allfold_dko_tameth", "subset1fold_mll3ko_tameth", "subset1fold_dko_tameth", "subset2fold_mll3ko_tameth", "subset2fold_dko_tameth"))

z$keytype <- factor(z$keytype, levels = c("allfold_mll3ko_tameth", "allfold_dko_tameth", "subset1fold_mll3ko_tameth", "subset1fold_dko_tameth"))
ggplot(data = z, aes(x=keytype, y=value, fill=type)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85") +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
 #ylim(c(-2,2)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- compare_means(value ~ keytype, data=z, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH")

#####
##### fold change scatterplot #####
genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X) # & fold_mll3ko_tameth > 0.5)
genesubset2 <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0.5)
genesubset3 <- df_lcpm_ert %>% filter(gene.name %in% genesubset$gene.name & gene.name %in% genesubset2$gene.name)

genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Wnt7b", "Cldn4", "Epcam")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

ggplot() +
  geom_point(data=genesubset, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth), color = "gray", size = 0.8) +
  geom_point(data=genesubset2, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth), color = "orange2", size = 0.8) +
  geom_text_repel(data = highlight, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  geom_smooth(data=genesubset, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth), method = lm, se=FALSE) + 
  xlim(-3,5) +
  ylim(-3,5) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X) # & fold_mll3ko_tameth > 0.5)
genesubset2 <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 > 0.5)
genesubset3 <- df_lcpm_ert %>% filter(gene.name %in% genesubset$gene.name & gene.name %in% genesubset2$gene.name)

genestohighlight <- c("Cldn6", "Tmem54", "Tacstd2", "Acsl5", "Wnt7b", "Cldn4", "Epcam")
highlight <- df_lcpm_ert %>% filter(gene.name %in% genestohighlight) #& (fold_mll3ko_tameth > 0.5 | fold_mll3ko_tameth < -0.5))

ggplot() +
  geom_point(data=genesubset, aes(x=fold_mll3ko_tameth16, y=fold_dko_tameth16), color = "gray", size = 0.4) +
  geom_point(data=genesubset2, aes(x=fold_mll3ko_tameth16, y=fold_dko_tameth16), color = "orange2", size = 0.4) +
  geom_text_repel(data = highlight, aes(x=fold_mll3ko_tameth16, y=fold_dko_tameth16, label=gene.name, size = 12), color = "black", force = 75, box.padding = 1)+
  geom_smooth(data=genesubset, aes(x=fold_mll3ko_tameth16, y=fold_dko_tameth16), method = lm, se=FALSE) + 
  xlim(-3,6) +
  ylim(-3,6) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X)
xyz <- as.data.frame(genesubset[,c("fold_mll3ko_tameth","fold_dko_tameth")])
x <- genesubset$fold_mll3ko_tameth
y <- genesubset$fold_dko_tameth
lm_eqn(filteredlcpm)
#y=0.32x - 0.16

genesubset <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X)
xyz <- as.data.frame(genesubset[,c("fold_mll3ko_tameth16","fold_dko_tameth16")])
x <- genesubset$fold_mll3ko_tameth16
y <- genesubset$fold_dko_tameth16
lm_eqn(filteredlcpm)
#y=0.45x+0.02

#####
##### refseq analysis #####

#UCSC reference processed and available in raw data folder
refseq <- read.table("~/Desktop/data/refseqanalysis/ucsc_refseq_TSS.txt")
colnames(refseq) <- c("transcript", "chr", "strand", "start", "end", "gene.symbol")

plus_strand <- refseq %>% filter(strand == "+")
plus_strand$tss_start <- plus_strand$start - 1000
plus_strand$tss_end <- plus_strand$start + 1000

minus_strand <- refseq %>% filter(strand == "-")
minus_strand$tss_start <- minus_strand$end - 1000
minus_strand$tss_end <- minus_strand$end + 1000

tss_regions <- rbind(minus_strand, plus_strand)
tss_regions2 <- as.data.frame(cbind(tss_regions$chr, tss_regions$tss_start, tss_regions$tss_end, tss_regions$gene.symbol))
colnames(tss_regions2) <- c("chr", "tss_start", "tss_end", "gene.symbol")
tss_regions2$tss_start <- sapply(tss_regions2$tss_start, as.numeric)
tss_regions2$tss_end <- sapply(tss_regions2$tss_end, as.numeric)
readcounts_tss <- read.table("~/Desktop/data_grhl/epigenomics/readcounts/readcounts_tss_ckoertall_1000.tab")

colnames(readcounts_tss) <- c("chr", "start", "end", 
                              "mll3ko_eth_k4m1", "mll3ko_tam_k4m1","dko_eth_k4m1", "dko_tam_k4m1", 
                              "mll3ko_eth_k4m2", "mll3ko_tam_k4m2","dko_eth_k4m2", "dko_tam_k4m2", 
                              "mll3ko_eth_k4m3", "mll3ko_tam_k4m3","dko_eth_k4m3", "dko_tam_k4m3",
                              "mll3ko_eth_k27a", "mll3ko_tam_k27a","dko_eth_k27a", "dko_tam_k27a", 
                              "mll3ko_eth_h33", "mll3ko_tam_h33","dko_eth_h33", "dko_tam_h33", 
                              "mll3ko_eth_rad21", "mll3ko_tam_rad21","dko_eth_rad21", "dko_tam_rad21",
                              "mll3ko_eth_ha", "mll3ko_tam_ha","dko_eth_ha", "dko_tam_ha", 
                              "mll3ko_eth_p300", "mll3ko_tam_p300","dko_eth_p300", "dko_tam_p300", 
                              "mll3ko_eth_igg", "mll3ko_tam_igg","dko_eth_igg", "dko_tam_igg")


readcounts_oi <- readcounts_tss
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:39] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)

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

readcounts_tss <- cbind(readcounts_oi2, fold_mll3ko_k4m1,fold_dko_k4m1,fold_mll3ko_k4m2,fold_dko_k4m2,
                        fold_mll3ko_k4m3, fold_dko_k4m3, fold_mll3ko_k27a, fold_dko_k27a, fold_mll3ko_h33, fold_dko_h33,
                        fold_mll3ko_rad21, fold_dko_rad21, fold_mll3ko_ha, fold_dko_ha,
                        fold_mll3ko_p300, fold_dko_p300, fold_mll3ko_igg, fold_dko_igg)

readcounts_tss2 <- left_join(readcounts_tss, tss_regions2, by = c("chr" = "chr", "start" = "tss_start", "end" = "tss_end"))
readcounts_tss3 <- distinct(readcounts_tss2, readcounts_tss2[36:38], .keep_all = TRUE)
readcounts_tss4 <- readcounts_tss3 %>% filter(gene.symbol %in% df_lcpm_ert$gene.name)

readcounts_tss3v2 <- readcounts_tss3 %>% replace(is.na(.), 0) %>% mutate(sum1 = rowSums(across("mll3ko_eth_k27a":"dko_tam_k27a")))
readcounts_tss3v3 <- readcounts_tss3v2 %>% group_by(gene.symbol) %>% dplyr::slice(which.max(sum1))
readcounts_tss4 <- readcounts_tss3v3 %>% filter(gene.symbol %in% df_lcpm_ert$gene.name)


#####
##### plot refseq tss counts #####

readcounts_tss_subset <- readcounts_tss4 %>% filter(gene.symbol %in% mll3kotameth8up$gene.name)

comparison <- c("fold_mll3ko_igg", "fold_dko_igg")
comparison <- c("fold_mll3ko_ha", "fold_dko_ha")
comparison <- c("fold_mll3ko_k4m1", "fold_dko_k4m1")
comparison <- c("fold_mll3ko_k4m2", "fold_dko_k4m2")
comparison <- c("fold_mll3ko_k4m3", "fold_dko_k4m3")
comparison <- c("fold_mll3ko_k27a", "fold_dko_k27a")
comparison <- c("fold_mll3ko_p300", "fold_dko_p300")

stats <- 0

comparison <- c("fold_mll3ko_igg", "fold_dko_igg")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q1 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-4,6)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH")

##
comparison <- c("fold_mll3ko_ha", "fold_dko_ha")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q2 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-4,6)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))
##
comparison <- c("fold_mll3ko_k4m1", "fold_dko_k4m1")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q3 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-2,4)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))
####
comparison <- c("fold_mll3ko_k4m2", "fold_dko_k4m2")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q4 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-2,4)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))
####
comparison <- c("fold_mll3ko_k4m3", "fold_dko_k4m3")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q5 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-2,4)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))
####
comparison <- c("fold_mll3ko_k27a", "fold_dko_k27a")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q6 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-2,4)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))
####
comparison <- c("fold_mll3ko_p300", "fold_dko_p300")
x <- gather(readcounts_tss4[,comparison])
x$kind <- "allactive"
x2 <- gather(readcounts_tss_subset[,comparison])
x2$kind <- "subset"
x <- rbind(x, x2)
x$key <- factor(x$key, levels = comparison)
x$keytype <- paste0(x$kind, x$key)
x$keytype <- factor(x$keytype, levels = c(unique(x$keytype)))
q7 <- ggplot(data = x, aes(x=keytype, y=value, fill=kind)) +
  scale_fill_manual(values=c("gray85", "orange2")) +
  geom_jitter(width=0.1, color = "gray85", size = 0.2) +
  stat_boxplot(geom = "errorbar", width = 0.2, color ="black") +
  geom_boxplot(width=0.9,outlier.shape=NA, color = "black") +
  ylim(c(-2,4)) +
  theme_bw() +
  theme(aspect.ratio = 2, text = element_text(size = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")
stats <- rbind(stats, compare_means(value ~ keytype, data=x, paired = FALSE, method = "wilcox.test", p.adjust.method = "BH"))


q1|q2|q3|q4|q5|q6|q7

q6

#####
##### venn diagram comparisons #####

mll3kowt8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt8hr$X & fold_mll3ko_wteth > 0.5) #23
mll3kotameth8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0.5) #101
dkotameth8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & fold_dko_tameth > 0.5) #9

nrow(mll3kowt8up %>% filter(gene.name %in% mll3kotameth8up$gene.name)) #0
nrow(mll3kowt8up %>% filter(gene.name %in% dkotameth8up$gene.name)) #0
nrow(mll3kotameth8up %>% filter(gene.name %in% dkotameth8up$gene.name)) #7
grid.newpage()                                        
draw.pairwise.venn(9,101,7)


mll3kowt16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt16hr$X & fold_mll3ko_wteth16 > 0.5) #4
mll3kotameth16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 > 0.5) #135
dkotameth16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & fold_dko_tameth16 > 0.5) #3
nrow(mll3kowt16up %>% filter(gene.name %in% mll3kotameth16up$gene.name)) #0
nrow(mll3kowt16up %>% filter(gene.name %in% dkotameth16up$gene.name)) #0
nrow(mll3kotameth16up %>% filter(gene.name %in% dkotameth16up$gene.name)) #7
grid.newpage()                                        
draw.pairwise.venn(3, 135, 2)

#relaxed gating
mll3kowt8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt8hr$X & fold_mll3ko_wteth > 0) #109
mll3kotameth8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0) #131
dkotameth8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko8hr$X & fold_dko_tameth > 0) #18

nrow(mll3kowt8up %>% filter(gene.name %in% mll3kotameth8up$gene.name)) #1
nrow(mll3kowt8up %>% filter(gene.name %in% dkotameth8up$gene.name)) #0
nrow(mll3kotameth8up %>% filter(gene.name %in% dkotameth8up$gene.name)) #11
grid.newpage()                                        
draw.pairwise.venn(18,131,11)

mll3kowt16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_ckowt16hr$X & fold_mll3ko_wteth16 > 0) #30
mll3kotameth16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko16hr$X & fold_mll3ko_tameth16 > 0) #445
dkotameth16up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_dko16hr$X & fold_dko_tameth16 > 0) #3
nrow(mll3kowt16up %>% filter(gene.name %in% mll3kotameth16up$gene.name)) #3
nrow(mll3kowt16up %>% filter(gene.name %in% dkotameth16up$gene.name)) #0
nrow(mll3kotameth16up %>% filter(gene.name %in% dkotameth16up$gene.name)) #3
grid.newpage()                                        
draw.pairwise.venn(3,445,3)

##cko8hr vs dox12hr

mll3kotameth8up <- df_lcpm_ert %>% filter(gene.name %in% diffgeneset_cko8hr$X & fold_mll3ko_tameth > 0.5)
diffgenesetdox12v0_up <- diffgeneset_doxsig12v0 %>% filter(foldchange12v0 > 1)
nrow(mll3kotameth8up %>% filter(gene.name %in% diffgeneset_doxsig12v0$gene.name )) #72
grid.newpage()                                        
draw.pairwise.venn(101,900,72)


####
#####
#####scatterplots at TSS #####

xx<- right_join(readcounts_tss4, mll3kotameth8up, by = c("gene.symbol" = "gene.name"), .keep_all)  
#filteredxx <- xx %>% filter(fold_mll3ko_tameth > 1 & fold_dko_tameth > 1)
filteredxx <- xx %>% filter(gene.symbol %in% dkotameth8up$gene.name & fold_dko_tameth > 0 & fold_mll3ko_tameth > 0)

e1 <- ggplot(data=xx, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_tameth, y=fold_dko_tameth), color = "orange2") +
  xlim(c(-2,5)) +
  ylim(c(-2,5)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e2 <- ggplot(data=xx, aes(x=fold_mll3ko_igg, y=fold_dko_igg)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_igg, y=fold_dko_igg), color = "orange2") +
  xlim(c(-5,10)) +
  ylim(c(-5,10)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")


e3 <- ggplot(data=xx, aes(x=fold_mll3ko_ha, y=fold_dko_ha)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_ha, y=fold_dko_ha), color = "orange2") +
  xlim(c(-2,10)) +
  ylim(c(-2,10)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e4 <- ggplot(data=xx, aes(x=fold_mll3ko_k4m1, y=fold_dko_k4m1)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_k4m1, y=fold_dko_k4m1), color = "orange2") +
  xlim(c(-1,3)) +
  ylim(c(-1,3)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e5 <- ggplot(data=xx, aes(x=fold_mll3ko_k4m2, y=fold_dko_k4m2)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_k4m2, y=fold_dko_k4m2), color = "orange2") +
  xlim(c(-1,1.5)) +
  ylim(c(-1,1.5)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e6 <- ggplot(data=xx, aes(x=fold_mll3ko_k4m3, y=fold_dko_k4m3)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_k4m3, y=fold_dko_k4m3), color = "orange2") +
  xlim(c(-1,1.5)) +
  ylim(c(-1,1.5)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e7 <- ggplot(data=xx, aes(x=fold_mll3ko_k27a, y=fold_dko_k27a)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_k27a, y=fold_dko_k27a), color = "orange2") +
  xlim(c(-2,3)) +
  ylim(c(-2,3)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

e8 <- ggplot(data=xx, aes(x=fold_mll3ko_p300, y=fold_dko_p300)) +
  geom_point() +
  geom_point(data=filteredxx, aes(x=fold_mll3ko_p300, y=fold_dko_p300), color = "orange2") +
  xlim(c(-1,3)) +
  ylim(c(-1,3)) +
  theme_bw() +
  theme(aspect.ratio = 1, text = element_text(linewidth = 10), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")


(e1|e2|e3|e4)/(e5|e6|e7|e8)


#####










