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


## Prep. workspace
if(!dir.exists("1.Ident")) dir.create("1.Ident")
setwd("1.Ident")

# Load cytoband
b_file <- "hg38_arm.rds"
if(!file.exists(b_file)){
	# Download
	zip_file <- tempfile(fileext = ".gz")
	download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz", zip_file, mode = "wb")
	band <- read.delim(zip_file, header = FALSE)
	colnames(band) <- c("chr", "start", "end", "name", "gieStain")
	band$chr <- gsub(x = band$chr, pattern = "chr0", replacement = "chr", fixed = TRUE)
	band$start <- band$start + 1
	band <- as(band, "GRanges")
	# Drop unmapped contigs
	seqlevels(band, pruning.mode = "coarse") <- paste0("chr", c(1:22, "X", "Y"))
	# Convert from cytobands to chr arms
	arm <- cytobandToArm(band)
	# Remove centromere
	arm <- arm[which(arm$type != "cen")]
	# Switch chromosome names
	GenomeInfoDb::seqlevelsStyle(arm) <- "Ensembl"
	saveRDS(arm, file = b_file)
}
b_file <- "hg38_arm.rds"
arm <- readRDS(file = b_file)


## Load the sc reference used in Redeconve 
sc_ref <- readRDS("../Redeconve_ref.rds")


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


	## ----------
	# Diagnostic plot: SME clusters
	p <- SpatialDimPlot(obj, group.by = "SME.cluster", stroke = 0)
	p <- p + scale_fill_manual(values = cluster_colors)
	ggsave(
		p, 
		width = 4, height = 3, 
		file = sprintf("%splot_cluster.png", s_name)
	)


	## ----------
	# Load CopyKAT predictions
	copykat <- readRDS(gsub(x = st_file, pat = ".rds", rep = "_copykat.rds"))
	# Attach CopyKAT predictions to the Seurat object
	obj$copykat <- with(
		copykat$prediction, 
		copykat.pred[match(colnames(obj), cell.names)]
	)
	if(anyNA(obj$copykat.pred))
		stop("Something went wrong")
	if(!is.factor(obj$copykat.pred))
		obj$copykat.pred <- factor(
			obj$copykat.pred, 
			levels = c("aneuploid", "diploid", "not.defined")
		)
	# Diagnostic plot: CopyKAT predictions
	p <- SpatialDimPlot(obj, group.by = "copykat.pred", stroke = 0)
	p <- p + scale_fill_manual(values = ploidy_colors)
	ggsave(
		p, 
		width = 4, height = 3, 
		file = sprintf("%splot_ploidy.png", s_name)
	)


	## ----------
	# Load cell fraction estimates
	ct_frac <- readRDS(gsub(x = st_file, pat = ".rds", rep = "_Redeconve.rds"))
	if(any(rownames(ct_frac) != sc_ref$barcodes))
		stop("Something went wrong")
	if(any(colnames(ct_frac) != colnames(obj)))
		stop("Something went wrong")

	# Aggregate cell fraction estimates: barcode level to lineage level
	ct_frac_split1 <- lapply(
		X = split(
			x = as.matrix(ct_frac), 
			f = sc_ref$annotations, 
			drop = TRUE
		), 
		FUN = matrix, 
		ncol = ncol(ct_frac)
	)
	ct_frac_merge1 <- do.call(
		what = rbind, 
		args = lapply(X = ct_frac_split1, FUN = colSums)
	)
	colnames(ct_frac_merge1) <- colnames(ct_frac)
	# Normalize cell fraction, so that they sum to 1 for each spot
	ct_frac_merge1 <- sweep(
		x = ct_frac_merge1, 
		MARGIN = 2, 
		STATS = colSums(ct_frac_merge1), 
		FUN = "/"
	)
	if(anyNA(ct_frac_merge1))
		stop("Something went wrong")

	# Aggregate cell fraction estimates: barcode level to sub-cluster level
	ct_frac_split2 <- lapply(
		X = split(
			x = as.matrix(ct_frac), 
			f = sc_ref$subcluster, 
			drop = TRUE
		), 
		FUN = matrix, 
		ncol = ncol(ct_frac)
	)
	ct_frac_merge2 <- do.call(
		what = rbind, 
		args = lapply(X = ct_frac_split2, FUN = colSums)
	)
	colnames(ct_frac_merge2) <- colnames(ct_frac)
	# Normalize cell fraction, so that they sum to 1 for each spot
	ct_frac_merge2 <- sweep(
		x = ct_frac_merge2, 
		MARGIN = 2, 
		STATS = colSums(ct_frac_merge2), 
		FUN = "/"
	)
	if(anyNA(ct_frac_merge2))
		stop("Something went wrong")
	# To avoid 0 * log2(0) = NaN in entropy estimation, set 0s to 1e-5
	ct_frac_merge2[ct_frac_merge2 < 1e-5] <- 1e-5
	if(!isTRUE(all.equal(
		ct_frac_merge1["Tumor", ], 
		ct_frac_merge2["LOXpos", ] + ct_frac_merge2["LOXneg", ], 
		tolerance = 1e-3
	))) stop("Something went wrong")
	saveRDS(as.data.frame(t(ct_frac_merge2)), file = sprintf("%scluster_subcluster.rds", s_name))

	# Diagnostic plot: major lineage fractions
	meta_cols <- sprintf("frac_%s", rownames(ct_frac_merge1))
	n_row <- length(meta_cols) %/% 3L + 1L
	obj <- AddMetaData(obj, t(ct_frac_merge1), col.name = meta_cols)
	p <- SpatialFeaturePlot(obj, features = meta_cols, stroke = 0, ncol = 3)
	ggsave(
		p, 
		width = 9, height = 3 * n_row, 
		file = sprintf("%splot_frac.png", s_name)
	)

	# Diagnostic plot: marker gene expression
	p <- SpatialFeaturePlot(
		obj, 
		features = c(
			"CA9", "PTPRC", "CD8A", "CD4", "CD79A", 
			"LYZ", "PECAM1", "ACTA2"
		), 
		stroke = 0, 
		ncol = 3
	)
	ggsave(
		p, 
		width = 9, height = 9, 
		file = sprintf("%splot_exp.png", s_name)
	)

	# Diagnostic plot: QC metrics
	p <- SpatialFeaturePlot(
		obj, 
		features = c("nCount_Spatial", "nFeature_Spatial"), 
		stroke = 0, 
		ncol = 2
	)
	ggsave(
		p, 
		width = 6, height = 3, 
		file = sprintf("%splot_QC.png", s_name)
	)


	## ----------
	# CAA
	# Convert values in copykat CNA matrix to copy number differences 
	# (in diploid)
	# Genomic positions to genomic regions
	gp <- GPos(
		seqnames = Rle(copykat$CNAmat$chrom), 
		pos = copykat$CNAmat$chrompos
	)
	gp <- sort(gp)
	gr <- lapply(
		X = split(x = gp, f = seqnames(gp)), 
		FUN = function(g){
			r <- as(g[-length(g)], "GRanges")
			end(r) <- start(g)[-1] - 1L
			return(r)
		}
	)
	gr <- do.call(c, unname(gr))
	copykat_res <- SummarizedExperiment(
		rowRanges = gr, 
		assays = list(
			diff = matrix(
				NA_real_, 
				nrow = length(gr), ncol = ncol(obj), 
				dimnames = list(NULL, colnames(obj))
			)
		)
	)
	CNAmat <- with(
		copykat, {
			tmp <- CNAmat
			colnames(tmp) <- prediction$cell.names[match(
				colnames(tmp), 
				rownames(prediction)
			)]
			data.matrix(tmp[, !is.na(colnames(tmp))]) * 2L
		}
	)
	if(anyNA(CNAmat))
		stop("Something went wrong")
	o <- findOverlaps(query = gp, subject = gr)
	if(length(unique(subjectHits(o))) != length(gr))
		stop("Something went wrong")
	assay(copykat_res)[subjectHits(o), colnames(CNAmat)] <- CNAmat[queryHits(o), ]
	# Call fraction of arms with CNA
	arm_frac <- regionAneuploidy(
		copykat_res, 
		region = arm, 
		assay = "diff", 
		baseline = 0, 
		ta = 0.1, 
		td = 0.1, 
		minoverlap = 0.7
	)
	arm_call <- regionCall(arm_frac, brlen = 0.5)
	if(!all(colnames(arm_call) == colnames(obj)))
		stop("Something went wrong")


	## ----------
	# Export plotting data
	to_plot <- cbind.data.frame(
		obj@meta.data[c(
			"orig.ident", 
			"nCount_Spatial", "nFeature_Spatial", 
			"SME.cluster"
		)], 
		obj@meta.data[meta_cols], 
		t(arm_call), 
		CAA = colSums(abs(arm_call), na.rm = TRUE)
	)
	saveRDS(to_plot, file = sprintf("%scluster_data.rds", s_name))

	# Clean
	rm(list = c("ct_frac"))
}


## Plot cell fraction and CNA
to_plot <- lapply(X = list.files(pattern = "cluster_data.rds"), FUN = readRDS)
to_plot <- do.call(rbind.data.frame, to_plot)
to_plot$orig.ident <- factor(to_plot$orig.ident)
# To reorder the colors of clusters
to_plot$SME.cluster <- factor(
	to_plot$SME.cluster, 
	levels = rev(sprintf(
		"niche%d", 
		seq_len(length(unique(to_plot$SME.cluster))) - 1
	))
)
cluster_colors <- setNames(
	cluster_colors_all[seq.int(nlevels(to_plot$SME.cluster))], 
	levels(to_plot$SME.cluster)
)

# Left
to_p1 <- data.frame(
	orig.ident = to_plot$orig.ident, 
	SME.cluster = to_plot$SME.cluster, 
	frac_Tumor = to_plot$frac_Tumor, 
	frac_immune = rowSums(data.matrix(to_plot[, sprintf(
		fmt = "frac_%s", 
		c("CD4T", "CD8T", "Macro", "Mono", "Neutrophil", "NK")
	)])), 
	frac_EC = to_plot$frac_EC, 
	frac_Fib = to_plot$frac_Fib
)
# frac(Tumor) + frac(Immu) + frac(EC) + frac(Fib) == 1
if(!isTRUE(all.equal(
	unname(rowSums(data.matrix(to_p1[, 3:6]))), 
	rep(1, nrow(to_p1))
))) stop("Something went wrong")
p1_list <- lapply(
	X = colnames(to_p1)[3:6], 
	FUN = function(i){
		p <- ggboxplot(
			data = to_p1, 
			x = "SME.cluster", y = i, 
			color = "SME.cluster", 
			outlier.shape = 20, 
			outlier.size = 1, 
			palette = cluster_colors, 
			orientation = "horizontal"
		)
		p <- p + facet_grid(orig.ident~., scales = "free", space = "free")
		p + labs(fill = NULL, x = NULL) + theme(legend.position = "none")
	}
)

# Middle
p2 <- ggviolin(
	data = to_plot, 
	x = "SME.cluster", 
	y = "nFeature_Spatial", 
	color = "SME.cluster", 
	outlier.shape = 20, 
	outlier.size = 1, 
	palette = cluster_colors, 
	orientation = "horizontal"
)
p2 <- p2 + facet_grid(orig.ident~., scales = "free", space = "free")
p2 <- p2 + labs(fill = NULL, x = NULL) + theme(legend.position = "none")

# Right: chr3
if(!all(to_plot[["chr3:p"]] %in% c(-1, 0, 1)))
	stop("Something went wrong")
p3 <- ggcrosstab(
	data = data.frame(
		SME.cluster = to_plot$SME.cluster, 
		orig.ident = to_plot$orig.ident, 
		chr3p = factor(
			c("Loss", "Neutral", "Gain")[to_plot[["chr3:p"]] + 2L], 
			levels = c("Gain", "Neutral", "Loss")
		)
	), 
	x = "SME.cluster", 
	fill = "chr3p", 
	orientation = "horizontal", 
	palette = c("Loss" = "blue", "Neutral" = "grey", "Gain" = "red"), 
	method = NULL
)
p3$layers[[2]] <- NULL
p3 <- p3 + facet_grid(orig.ident~., scales = "free", space = "free")
p3 <- p3 + labs(y = "Pct. spots with 3p loss") + theme(legend.position = "bottom")
# Right: chr5
if(!all(to_plot[["chr5:q"]] %in% c(-1, 0, 1)))
	stop("Something went wrong")
p4 <- ggcrosstab(
	data = data.frame(
		SME.cluster = to_plot$SME.cluster, 
		orig.ident = to_plot$orig.ident, 
		chr5q = factor(
			c("Loss", "Neutral", "Gain")[to_plot[["chr5:q"]] + 2L], 
			levels = c("Gain", "Neutral", "Loss")
		)
	), 
	x = "SME.cluster", 
	fill = "chr5q", 
	orientation = "horizontal", 
	palette = c("Loss" = "blue", "Neutral" = "grey", "Gain" = "red"), 
	method = NULL
)
p4$layers[[2]] <- NULL
p4 <- p4 + facet_grid(orig.ident~., scales = "free", space = "free")
p4 <- p4 + labs(y = "Pct. spots with 5q gain") + theme(legend.position = "bottom")

# Merge plots
p <- wrap_plots(c(p1_list, list(p2, p3, p4)), nrow = 1, guides = "collect")
ggsave(p, file = "Niche_statistics_plot.png", width = 12, height = 6, bg = "white")
