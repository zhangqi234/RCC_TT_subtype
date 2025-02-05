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
if(!dir.exists("3.Scoring")) dir.create("3.Scoring")
setwd("3.Scoring")


## Partition spots into discrete regions based on tumor niche scores
parts <- delimitGeoRegions(d = 2L, size = 5L)
# Save and plot
saveRDS(parts[c("breaks", "levels")], file = "partition_model.rds")
ggsave(parts, file = "partition_model.png", width = 8, height = 3.2, bg = "white")
ggsave(parts, file = "partition_model.pdf", width = 8, height = 3.2, bg = "transparent")

# Plot the color scale
p_scale <- ggplot() + geom_raster(
	data = data.frame(
		x = seq.int(length(parts$levels)), 
		y = 1L, 
		z = sort(parts$levels, decreasing = FALSE)
	), 
	mapping = aes(x = x, y = y, fill = z), 
	interpolate = FALSE, 
	position = "identity"
) + viridis::scale_fill_viridis(
	option = "F", 
	direction = -1, 
	limits = c(-0.01, 1.01)
)
p_scale <- p_scale + scale_x_continuous(
	name = "Relative distance to tumor core", 
	expand = expansion(mult = c(0, 0)), 
	breaks = p_scale$layers[[1]]$data$x - 0, 
	labels = {
		l <- character(length = nrow(p_scale$layers[[1]]$data))
		idx <- match(parts$breaks, p_scale$layers[[1]]$data$z)
		l[idx] <- names(parts$breaks)
		l
	}
) + scale_y_continuous(
	expand = expansion(add = c(0, 0))
) + theme(
	legend.position = "none"
)


## Load tumor core definition
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
	# Plot tumor core gradient
	obj$tc <- ifelse(as.character(obj$SME.cluster) %in% tc[[s_name]], "core", "peri")
	obj$frac_LOXpos <- meta$LOXpos
	obj$frac_Fib5 <- meta$F05_Fib_COL1A1
	obj$frac_Macro5 <- meta$M05_Macro_APOC1
	# Partition spots based on relative distance to tumor core niches
	obj$region_score2 <- getGeoScores(object = obj, idents = (obj$tc == "core"), d = 2L, compare = FALSE)
	if(anyNA(obj$region_score2)) stop("Something went wrong")
	obj$region_score2 <- round(obj$region_score2, digits = 3L)
	if(any(obj$region_score2 < 0) || any(obj$region_score2 > 1)) stop("Something went wrong")
	obj$region_partition2 <- cut(obj$region_score2, breaks = c(-Inf, parts$levels), right = TRUE, include.lowest = TRUE)
	if(anyNA(obj$region_partition2)) stop("Something went wrong")
	tmp <- table(obj$region_partition2)
	idx <- which(tmp < 10L)
	if(length(idx) != 0L) warning("<10 spots with TC scores ", paste(names(tmp)[idx], collapse = ", "))
	# Color spots by TC region score
	p <- SpatialDimPlot(obj, group.by = "region_score2", pt.size = 1, pt.shape = "hexagon", greyscale = TRUE)
	p$data$region_score2 <- as.numeric(as.character(p$data$region_score2))
	p <- p + viridis::scale_fill_viridis(
		option = "F", 
		direction = -1, 
		limits = c(-0.01, 1.01), 
		breaks = parts$breaks, 
		labels = names(parts$breaks), 
		name = "Relative distance to tumor core"
	)
	p$layers[[2]]$aes_params$alpha <- NULL
	p$mapping$alpha <- quote(region_score2)
	p <- p + scale_alpha_continuous(range = c(0.6, 0.9), trans = "reverse", guide = "none")
	# Overlay (LOX+ tumor, THBS2+ fib and GPNMB+ macro) or, just CD8+ T cells?
	CD8_exp <- FetchData(obj, c("CD8A", "CD8B"), cells = rownames(p$data), layer = "data")
	p <- p + geom_point(
		data = p$data[which(rowSums(CD8_exp) != 0), ], 
		inherit.aes = FALSE, 
		mapping = aes(x = imagecol, y = imagerow), 
		colour = "#00DDFF", 
		fill = NA, 
		shape = 16, 
		alpha = 1, 
		size = p$layers[[2]]$aes_params$size
	)
	ggsave(p, file = sprintf("%splot_partition.png", s_name), width = 5.5, height = 2.5, bg = "white")
	ggsave(p, file = sprintf("%splot_partition.pdf", s_name), width = 5.5, height = 2.5, bg = "transparent")


	## ----------
	# Estimate co-localization scores
	# Get TIB co-localization scores (as weighted evenness) for each spot
	ct_tmp <- data.matrix(data.matrix(meta[, c("LOXpos", "F05_Fib_COL1A1", "M05_Macro_APOC1")]))
	ct_tmp2 <- ct_tmp / rowSums(ct_tmp)
	if(anyNA(ct_tmp)) stop("Something went wrong")
	obj$TIB_score2 <- -rowSums(ct_tmp2 * log2(ct_tmp2)) / log2(3) * rowSums(ct_tmp)
	if(anyNA(obj$TIB_score2)) stop("Something went wrong")

	# Diagnostic plot: co-localization scores
	p <- SpatialFeaturePlot(obj, features = "TIB_score2", stroke = 0)
	p <- p + theme(legend.position = "right")
	ggsave(p, file = sprintf("%splot_TIB2.png", s_name), width = 4, height = 3, bg = "white")

	# Plot co-localization scores and pct. CD8+ spots at each partition
	to_plot <- obj@meta.data[, c("region_partition2", "TIB_score2")]
	p <- ggline(
		data = to_plot, 
		x = "region_partition2", 
		y = "TIB_score2", 
		color = "red", 
		size = 0.2, 
		linewidth = 0.2, 
		point.size = 0.2, 
		add = "mean_se", 
		add.params = list(width = 0.3), 
		legend = "right"
	)
	CD8_exp <- FetchData(obj, c("CD8A", "CD8B"), cells = colnames(obj), layer = "data")
	to_plot2 <- aggregate(
		x = data.frame(CD8_exp_frac = rowSums(CD8_exp) != 0), 
		by = data.frame(region_partition2 = obj$region_partition2), 
		FUN = mean
	)
	CD8_scale <- 3
	if(max(to_plot2$CD8_exp_frac) * CD8_scale > 1)
		stop("CD8_scale is too large")
	p2 <- ggline(
		data = to_plot2, 
		x = "region_partition2", 
		y = "CD8_exp_frac", 
		color = "#00DDFF", 
		size = 0.2, 
		linewidth = 0.2, 
		point.size = 0.2, 
		legend = "right"
	)
	p2$mapping$y <- rlang::quo_set_expr(p2$mapping$y, expr = quote(CD8_exp_frac * CD8_scale))
	for(i in seq_along(p2$layers)) p2$layers[[i]]$data <- p2$data
	for(i in seq_along(p2$layers)) p2$layers[[i]]$mapping <- structure(c(p2$layers[[i]]$mapping, p2$mapping), class = "uneval")
	p$layers <- c(p$layers, p2$layers)
	p <- p + scale_x_discrete(
		expand = expansion(add = c(0.5, 0.5))
	) + scale_y_continuous(
		limits = c(0, 1), 
		expand = expansion(add = c(0, 0)), 
		sec.axis = sec_axis(trans = ~ ./CD8_scale, name = "Pct. CD8+ spots")
	) + labs(
		y = "Avg. co-localization strength"
	)
	ggsave(p, file = sprintf("%splot_trend_diagnostics.png", s_name), width = 15, height = 5, bg = "white")
	for(i in seq_along(p$layers)) if(is(p$layers[[i]]$geom, "GeomPoint"))
		p$layers[[i]]$aes_params$size <- 0.5
	p <- p + theme(axis.title.y.left = element_text(colour = "red"), axis.title.y.right = element_text(colour = "#00DDFF"))
	p_final <- wrap_plots(p, p_scale, ncol = 1, heights = c(10, 1)) & theme(plot.margin = unit(rep(0, 4), "pt"))
	ggsave(p_final, file = sprintf("%splot_trend_final.png", s_name), width = 2.5, height = 2, bg = "white")
	ggsave(p_final, file = sprintf("%splot_trend_final.pdf", s_name), width = 2.5, height = 2, bg = "transparent")

	## ----------
	# Compare pct. CD8+ spots: tumor nest v.s. peri-nest
	dat <- data.frame(
		tc = ifelse(obj$region_score2 <= parts$breaks["front"], yes = "core", no = "peri"), 
		CD8_exp = rowSums(FetchData(obj, c("CD8A", "CD8B"), cells = colnames(obj), layer = "data"))
	)
	dat$tc <- factor(dat$tc, levels = c("peri", "core"))
	dat$CD8_exp <- dat$CD8_exp != 0
	# Odds ratio and 95% confidence interval
	y <- glm(CD8_exp ~ tc, data = dat, family = binomial())
	or <- exp(coef(y)["tccore"])
	conf <- exp(confint.default(y)["tccore", ])
	pval <- summary(y)$coefficients["tccore", "Pr(>|z|)"]
	p <- forestmodel::forest_model(y)
	ggsave(p, file = sprintf("%splot_OR.png", s_name), width = 5, height = 0.8, bg = "white")
	ggsave(p, file = sprintf("%splot_OR.pdf", s_name), width = 5, height = 0.8, bg = "transparent")
}
