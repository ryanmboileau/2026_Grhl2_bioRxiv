library(tidyverse)
library(ggplot2)
library(ggpubr)

##import deeptools matrix output with firstline/header already removed. remove regions, isolate bin counts
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_testct16igg/grhl2_formsites_nf_ert.1070607', sep = '\t')
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_testct14ha/grhl2_formsites_nf_ert.1070600', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_testct14k4me1/grhl2_formsites_nf_ert.1070607', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_testct14k4me2/grhl2_formsites_nf_ert.1070608', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_testct14k27a/grhl2_formsites_nf_ert.1070607', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ert/new_GRHL2sites_nf_p300/grhl2_formsites_nf_ert.1070607', sep = '\t', skip = 1)
# deeptools_matrix <- tail(deeptools_matrix, 332)
# deeptools_matrix2 <- deeptools_matrix[1:6]
# df_bincounts_raw <- deeptools_matrix[7:ncol(deeptools_matrix)]

# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_testct16igg/grhl2_formsites_nf_ectopicert.1100744', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_testct14ha/grhl2_formsites_nf_ectopicert.1100738', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_testct14k4me1/grhl2_formsites_nf_ectopicert.1100744', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_testct14k4me2/grhl2_formsites_nf_ectopicert.1100756', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_testct14k27a/grhl2_formsites_nf_ectopicert.1100744', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf_ectopicert/new_GRHL2sites_nf_p300/grhl2_formsites_nf_ectopicert.1100744', sep = '\t', skip = 1)
# deeptools_matrix2 <- deeptools_matrix[1:6]
# df_bincounts_raw <- deeptools_matrix[7:ncol(deeptools_matrix)]


# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_grhl2/grhl2_formsites_nf.1140817', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_k4m1/grhl2_formsites_nf.1140817', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_k27a/grhl2_formsites_nf.1140811', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_ct17k27a/grhl2_formsites_nf.1140811', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_ct17igg/grhl2_formsites_nf.1140854', sep = '\t', skip = 1)
# deeptools_matrix <- read.csv('~/Desktop/data_grhl/epigenomics/heatmaps2/grhl2_formsites_nf/new_GRHL2sites_nf_p300/grhl2_formsites_nf.1140854', sep = '\t', skip = 1)

deeptools_matrix <- tail(deeptools_matrix, 332)
deeptools_matrix2 <- deeptools_matrix[1:6]
df_bincounts_raw <- deeptools_matrix[7:ncol(deeptools_matrix)]



##specify number of samples, binsize and range used in deeptools output
samples <- 4
bins <- 40
range <- 3000 

#specify number of total bins to use for statistic test, need even number
binsforstats <- 2

##determine what column numbers to compare with 0 indexing
sample_number <- c(1,3)


raw_array=list()
mean_array=list()
ci_high_array=list()
ci_low_array=list()

for(i in sample_number) {
  #decompose large input matrix into individual samples and their respective bins, calculate column average down each bin
  #to prepare for moving window averaging
  sample_n <- 2 * range / bins 
  bin_n <- (-1 * range) + (bins * seq(sample_n))
  
  leftcol <- (sample_n*i) + 1
  rightcol <- leftcol + sample_n - 1
  df_bincounts_targetn <- df_bincounts_raw[leftcol:rightcol]
  colnames(df_bincounts_targetn) <- bin_n
  df_bincounts_targetn_mean <- colMeans(df_bincounts_targetn[sapply(df_bincounts_targetn, is.numeric)])
  

  raw_array[[paste0(i)]] <- as.data.frame(df_bincounts_targetn)

  #conduct moving window averaging using n total bins
  moving_average <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
  df_bincounts_targetn_mean_smooth <- as.vector(moving_average(df_bincounts_targetn_mean))
  
  df_bincounts_targetn_mean_smooth2 <- as.data.frame(df_bincounts_targetn_mean_smooth)
  df_bincounts_targetn_mean_smooth2$bin <- bin_n
  mean_array[[paste0(i)]] <- as.data.frame(df_bincounts_targetn_mean_smooth2)
  
  #calculate confidence interval threshold for raw inputs at each bin, create plus/minus conf int values
  confint_threshold <- 0.95
  df_bincounts_targetn_sd <- sapply(df_bincounts_targetn,function(x) (sd(x)*(confint_threshold))/sqrt(length(x)))
  
  df_bincounts_targetn_sd_high <- df_bincounts_targetn_mean + df_bincounts_targetn_sd
  df_bincounts_targetn_sd_high2 <- as.data.frame(df_bincounts_targetn_sd_high)
  df_bincounts_targetn_sd_high2$bin <- c(seq(1,nrow(df_bincounts_targetn_sd_high2)))
  ci_high_array[[paste0(i)]] <- as.data.frame(df_bincounts_targetn_sd_high2)
  
  df_bincounts_targetn_sd_low <- df_bincounts_targetn_mean - df_bincounts_targetn_sd
  df_bincounts_targetn_sd_low2 <- as.data.frame(df_bincounts_targetn_sd_low)
  df_bincounts_targetn_sd_low2$bin <- c(seq(1,nrow(df_bincounts_targetn_sd_low2)))
  ci_low_array[[paste0(i)]] <- as.data.frame(df_bincounts_targetn_sd_low2)
  
  
}

ggplot() +
  geom_line(data=mean_array[[1]], aes(y=df_bincounts_targetn_mean_smooth, x=bin), size=2) +
  geom_line(data=mean_array[[2]], aes(y=df_bincounts_targetn_mean_smooth, x=bin), size=2, color = "dodgerblue2") +
  theme_classic() +
  theme(aspect.ratio = 1)



#####
#####comparison statistics#####

y <- bind_cols(mean_array, .id="column_label")
x <- bind_rows(mean_array, .id = "column_label")
x <- bind_rows(raw_array, .id = "column_label")

z <- binsforstats / 2 * bins
xx <- x %>% select(column_label ,toString(-z):toString(z))
xx <- xx %>% mutate(mean_sel = rowMeans(select(.,toString(-z):toString(z))))
xx$mean_sel <- as.numeric(xx$mean_sel)
xx$bin <- as.numeric(0)


compare_means(mean_sel ~ column_label, data = xx, paired = TRUE, method = "wilcox.test")
stats<-compare_means(mean_sel ~ column_label, data = xx, paired = TRUE, method = "wilcox.test")


#####