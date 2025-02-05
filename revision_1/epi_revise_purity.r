#!/usr/bin/Rscript
# Author: Tan Yezhen & Zhang Qi
# Date: 2024/10/29


library(ggplot2)
library(ggpubr)
library(SummarizedExperiment)
TT_color <- c("TT1" = "#ff7f0f", "TT2" = "#3cb7cc")
dataset_color <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[seq(from = 3, to = 7, length.out = 3)]


## Plotting
obj_called_PUTH <- readRDS(file = "../CNA/4.plot.called.rds")
obj_called_EGA <- readRDS(file = "../CNA_EGA/4.plot.called.rds")
pdata_PUTH <- as.data.frame(colData(obj_called_PUTH))
pdata_EGA <- as.data.frame(colData(obj_called_EGA))
stopifnot(identical(colnames(pdata_PUTH), colnames(pdata_EGA)))
pdata_PUTH$dataset <- "PUTH"
pdata_EGA$dataset <- "EGA"
combined_pdata <- rbind.data.frame(pdata_PUTH, pdata_EGA)
combined_pdata$group <- paste(combined_pdata$Cluster, combined_pdata$Tumour_type, sep = "_")
# Samples to use are those have ASCAT profiles and cluster labels
s_use <- rownames(combined_pdata)[combined_pdata$Include & (!is.na(combined_pdata$Cluster))]

# Plot sample-level data
#for(n in c("purity", "ploidy", "goodnessOfFit", "LOH", "GI", "CAA", "wGII")){
n <- "purity"
p <- ggbarplot(
	data = combined_pdata[s_use, c("Patient", "dataset", "group", n)], 
	x = "group", y = n, color = "black", fill = "transparent", 
	add = c("mean_se"), add.params = list(color = "black", width = 0.3),
	legend = "right", width = 0.7, orientation = "horizontal")
p <- p + geom_jitter(
	data = combined_pdata[s_use, c("Patient", "dataset", "group", n)], 
	mapping = aes(x = group, y = .data[[n]], color = dataset), 
	stroke = 0, shape = 16, size = 1.6, width = 0.25, height = 0, 
	show.legend = TRUE)
p <- p + stat_compare_means(
	comparisons = list(c("TT1_TT", "TT2_TT"), c("TT1_PT", "TT2_PT")), 
	#step.increase = 0.09, tip.length = 0.03, vjust = -0.015, 
	label = "p.format", 
	method = "wilcox.test")
p <- p + scale_color_manual(values = dataset_color[c(1, 3)])
p <- p + scale_x_discrete(limits = c("TT2_PT", "TT1_PT", "TT2_TT", "TT1_TT"))
p <- p + theme(panel.border = element_blank(), aspect.ratio = 0.75)
p$layers <- p$layers[c(3, 1, 2, 4)]

min_y <- 0.15
if(min(p$data[[n]], na.rm = TRUE) < min_y) stop("Something went wrong", rep("!\n", 10))
max_y <- 1
if(max(p$data[[n]], na.rm = TRUE) > max_y) stop("Something went wrong", rep("!\n", 10))
p <- p + coord_flip(ylim = c(min_y, max_y), clip = "on")
p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 1), breaks = seq(0, 1, by = 0.2), position = "right")

ggsave(p, file = "epi_revise_purity_WES.png", width = 4, height = 2, bg = "white")
ggsave(p, file = "epi_revise_purity_WES.pdf", width = 4, height = 2, bg = "transparent")
#}

