#!/usr/bin/Rscript
# Author: Tan Yezhen & Zhang Qi
# Date: 2024/02/02

library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(patchwork)
library(Seurat)
library(SummarizedExperiment)
source("../helper_functions/SpatialRegions.R")
source("../helper_functions/SpatialBorders.R")



## Constants
# Color definitions
ploidy_colors <- c(aneuploid = "#E41A1C", diploid = "#377EB8", not.defined = "grey")
LOX_colors <- c(pos = "#9200af", neg = "#6787f3")
cluster_colors_all <- RColorBrewer::brewer.pal(name = "Set3", n = 12)[-9]

# Cell type definitions
ct_list <- c(
	"B", "CD4T", "CD8T", "DC", "EC", 
	"Fib_other", "F05_Fib_COL1A1", "Tumor", "Macro_other", 
	"M05_Macro_APOC1", "Mast", "Mono", "NK", "Neutrophil", 
	"Plasma"
)
ct_TIB <- c(
	"F05_Fib_COL1A1", "M05_Macro_APOC1"
)
ct_TIB2 <- c(
	"F05_Fib_COL1A1", "M05_Macro_APOC1", "LOXpos"
)


## Create working space
if(!dir.exists("2.Plot")) dir.create("2.Plot")
setwd("2.Plot")


## Define tumor core
# Average CAA, 3p-, 5q+ in each niche & tissue section
sapply(
	simplify = FALSE, 
	X = list.files(pattern = "cluster_data.rds", path = "../1.Ident", full.names = TRUE), 
	FUN = function(f){
		tab <- readRDS(f)
		tab_summary <- lapply(
			X = split(
				x = factor(c("Loss", "Neutral", "Gain")[tab[["chr3:p"]] + 2L], levels = c("Gain", "Neutral", "Loss")), 
				f = tab$SME.cluster
			), 
			FUN = table
		)
		tab_summary <- do.call(rbind, tab_summary)
		round(tab_summary / rowSums(tab_summary), digits = 3) * 100
	}
)
# Classify a niche as a tumor nest:
# ((high 3p loss) and (high tumor purity)) or (with necrosis)
# $p2620_cluster_data.rds
#        Gain Neutral Loss
# niche0  0.1    74.6 25.2
# niche1  1.3    95.1  3.6
# niche2  2.3    88.7  9.0
# niche3  0.6    61.9 37.5
# niche4 16.0    84.0  0.0
# $p2806_cluster_data.rds
#        Gain Neutral Loss
# niche0  1.1    98.3  0.6
# niche1  3.1    85.2 11.7
# niche2  0.6    94.4  5.0
# niche3  0.0    63.2 36.8
# niche4  7.7    84.6  7.7
# niche5  5.0    91.2  3.8
# $p4828_cluster_data.rds
#        Gain Neutral Loss
# niche0    0   100.0  0.0
# niche1    0   100.0  0.0
# niche2    0    98.2  1.8
# niche3    0    90.2  9.8
# niche4    0    87.4 12.6
# niche5    0    99.4  0.6
# niche6    0   100.0  0.0
# niche7    0    97.0  3.0
# And make this 'DefTumorCore.txt' file accordingly.
# tc <- list(section_A = c("niche0", "niche5"), ...)
source("../DefTumorCore.txt")


## Load and process st data
st_list <- c(
	# Some Seurat Visium rds files here
)
for(st_file in st_list){
	# Load st data
	obj <- readRDS(st_file)
	DefaultAssay(obj) <- "SCT"
	# Get sample name
	s_name <- obj@images[[1L]]@key
	# Load previously saved plotting data
	meta <- readRDS(sprintf("../1.Ident/%scluster_data.rds", s_name))
	if(!all(rownames(meta) == colnames(obj))) stop("Something went wrong")


	## ----------
	# Plot H&E and tumor core boundary
	obj$tc <- ifelse(as.character(obj$SME.cluster) %in% tc[[s_name]], "core", "peri")
	if(anyNA(obj$tc)) stop("Something went wrong")
	p_added <- addSpatialBorders(
		Seurat::SpatialDimPlot(obj, group.by = "tc"), 
		highlights = list("core"), 
		idents = obj$tc, 
		size = 0.2, 
		shape = 16
	)
	p <- SpatialDimPlot(obj, group.by = "tc")
	p$layers[[2]] <- p_added$layers[[2]]
	ggsave(p, file = sprintf("%splot_core.png", s_name), width = 2.5, height = 2.5, bg = "white")
	ggsave(p, file = sprintf("%splot_core.pdf", s_name), width = 2.5, height = 2.5, bg = "transparent")


	## ----------
	# Plot epi frac and 3p- frac at niche level
	tmp <- meta[, grepl(x = colnames(meta), pattern = "^frac_", perl = TRUE) & sapply(meta, is.numeric)]
	to_plot1 <- structure(
		lapply(
			X = tmp, 
			FUN = function(j) unsplit(
				f = meta$SME.cluster, 
				lapply(split(x = j, f = meta$SME.cluster), FUN = mean)
			)
		), 
		class = "data.frame", 
		row.names = rownames(meta)
	)
	obj$frac <- to_plot1$frac_Tumor
	p <- SpatialDimPlot(obj, group.by = "frac", pt.size = 1, pt.shape = "hexagon", greyscale = TRUE)
	p$data$frac <- as.numeric(as.character(p$data$frac))
	p$data$SME.cluster <- meta[rownames(p$data), "SME.cluster"]
	p$scales$scales[sapply(p$scales$scales, function(s) any(s$aesthetics == "fill"))] <- NULL
	p <- p + scale_fill_gradientn(colors = RColorBrewer::brewer.pal(name = "YlOrRd", n = 9), limits = c(0, 1))
	#
	tmp <- meta[, "chr3:p", drop = TRUE] == -1
	to_plot2 <- structure(
		unsplit(
			f = meta$SME.cluster, 
			lapply(split(x = tmp, f = meta$SME.cluster), FUN = mean)
		), 
		names = rownames(meta)
	)
	obj$chr3p_loss <- to_plot2
	pos_df <- SpatialDimPlot(obj, group.by = "chr3p_loss", label = TRUE)$layers[[3]]$data
	pos_df$chr3p_loss <- as.numeric(as.character(pos_df$chr3p_loss))
	pos_df$chr3p_g_or_n <- 1 - pos_df$chr3p_loss
	p <- p + ggnewscale::new_scale_fill()
	p <- p + scatterpie::geom_scatterpie(
		data = pos_df, 
		mapping = aes(x = imagecol, y = imagerow, group = color), 
		color = "black", 
		cols = c("chr3p_loss", "chr3p_g_or_n"), 
		pie_scale = 3, 
		size = 0.3
	)
	p <- p + scale_fill_discrete(type = c("chr3p_loss" = "black", "chr3p_g_or_n" = "white"))
	p <- p + labs(fill_new = "Pct. tumor cell", fill = "Pct. 3p loss")
	ggsave(p, file = sprintf("%splot_frac_avg.png", s_name), width = 4, height = 2.5, bg = "white")
	ggsave(p, file = sprintf("%splot_frac_avg.pdf", s_name), width = 4, height = 2.5, bg = "transparent")
}


