##Figure 1 GRHL2-ERT proof of principle##
library(ggpubr)
library(patchwork)

#####
#####ert controls#####
#link to file in raw data folder
tableoi <- read.delim("~/Desktop/Kevinswork/ert_controls8hr_R.txt", header=T, sep="\t", row.names = NULL, check.name=FALSE,
                 stringsAsFactor=FALSE)

tableoi$samplecondition <- paste0(tableoi$sample, tableoi$condition)
tableoi$samplecondition <- factor(tableoi$samplecondition, levels = c("v6etoh","v6tam","ertetoh", "erttam"))

a1 <- ggplot(data=tableoi, aes(x=samplecondition, y=grhl2, color=condition)) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  stat_summary(fun=mean, geom="point", color="black") +
  geom_point() +
  ylab("") +
  xlab("") + geom_hline(yintercept=0) +
  scale_color_manual(values=c("gray30", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")

a2 <- ggplot(data=tableoi, aes(x=samplecondition, y=cldn6, color=condition)) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  stat_summary(fun=mean, geom="point", color="black") +
  geom_point() +
  ylab("") +
  xlab("") + geom_hline(yintercept=0) +
  scale_color_manual(values=c("gray30", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")

a1 + a2

######
######ert timecourse#####

#link to file in raw data folder
tableoi <- read.delim("~/Desktop/Kevinswork/ert_timecourse_R.txt", header=T, sep="\t", row.names = NULL, check.name=FALSE,
                      stringsAsFactor=FALSE)

tableoi_exon <- tableoi %>% filter(type == "exon_exon") 

b1 <- ggplot(data=tableoi_exon, aes(x=time, y=(cldn6))) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  stat_summary(fun=mean, geom="point", color="black") +
  geom_point() +
  ylab("") +
  xlab("") + 
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks = c(tableoi_exon$time)) +
  stat_summary(fun=mean, geom="line", colour="black") +
  #scale_color_manual(values=c("gray", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 0.5, text = element_text(size = 20), legend.position = "none")

b2 <- ggplot(data=tableoi_exon, aes(x=time, y=(wnt7))) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  stat_summary(fun=mean, geom="point", color="black") +
  geom_point() +
  ylab("") +
  xlab("") + 
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks = c(tableoi_exon$time)) +
  stat_summary(fun=mean, geom="line", colour="black") +
  #scale_color_manual(values=c("gray", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 0.5, text = element_text(size = 20), legend.position = "none")

##
tableoi_intron <- tableoi %>% filter(type == "intron_exon") 

b3 <- ggplot(data=tableoi_intron, aes(x=time, y=(cldn6))) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  #stat_summary(fun=mean, geom="point", color="black") +
  geom_point(shape = 1) +
  ylab("") +
  xlab("") + 
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks = c(tableoi_intron$time)) +
  stat_summary(fun=mean, geom="line", colour="black", linetype = "dashed") +
  #scale_color_manual(values=c("gray", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 0.5, text = element_text(size = 20), legend.position = "none") +
  scale_y_continuous(position = "right")

b4 <- ggplot(data=tableoi_intron, aes(x=time, y=(wnt7))) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  #stat_summary(fun=mean, geom="point", color="black") +
  geom_point(shape = 1) +
  ylab("") +
  xlab("") + 
  geom_hline(yintercept=0) +
  scale_x_continuous(breaks = c(tableoi_intron$time)) +
  stat_summary(fun=mean, geom="line", colour="black", linetype = "dashed") +
  #scale_color_manual(values=c("gray", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 0.5, text = element_text(size = 20), legend.position = "none") +
  scale_y_continuous(position = "right")

b1 + b2 + b3 + b4


######
######ert chipqpcr ######
#link to file in raw data folder
tableoi <- read.delim("~/Desktop/data_grhl/ert_chipqpcr_Rv2.txt", header=T, sep="\t", row.names = NULL, check.name=FALSE,
                      stringsAsFactor=FALSE)

gather_tableloi <- gather(data=tableoi,key="key", value="value", 3:7)
gather_tableloi$samplekey <- paste0(gather_tableloi$sample, gather_tableloi$key)
gather_tableloi$samplekey <- factor(gather_tableloi$samplekey, levels = c(unique(paste0(gather_tableloi$sample, gather_tableloi$key))))

ggplot(data=gather_tableloi, aes(x=samplekey, y=value)) +
#ggplot(data=gather_tableloi, aes(x=key, y=value, group = key, color = sample)) +
  geom_point(color = "grey") +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=0.5) +
  stat_summary(fun=mean, geom="point", color="black") +
  ylab("") +
  xlab("") + 
  geom_hline(yintercept=0) +
  #scale_color_manual(values=c("gray30", "darkorange1")) +
  theme_classic() +
  theme(aspect.ratio = 0.5, text = element_text(size = 8), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

stats <- compare_means(value ~ sample, data = gather_tableloi, method="t.test", paired=TRUE, group.by = c("key"), p.adjust.method="BH")

gather_tableloi2 <- gather_tableloi %>% filter(key == "nanog_prom_E3")
gather_tableloi2 <- gather_tableloi %>% filter(key == "cdh1_2")
gather_tableloi2 <- gather_tableloi %>% filter(key == "dsp_2")
gather_tableloi2 <- gather_tableloi %>% filter(key == "cldn6")
gather_tableloi2 <- gather_tableloi %>% filter(key == "wnt7b")
stats <- compare_means(value ~ sample, data = gather_tableloi2, method="t.test", paired=TRUE, group.by = c("key"), p.adjust.method="BH")




######






