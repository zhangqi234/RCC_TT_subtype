#!/usr/bin/Rscript
# Author: Tan Yezhen & Zhang Qi
# Date: 2024/11/03

options(future.globals.maxSize=160*1024^3)
library(monocle)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(SingleCellExperiment)
library(Seurat)

# Tumor + renal epithelial cells
DTTEpi <- readRDS("")
# Renal epithelial cells, sub-clustered into proximal tubules, intercalated 
# cells ...
DTANTEpi <- readRDS("")

# Remove normal renal epithelial cells other than proximal tubules.
del_barcode <- colnames(DTANTEpi)[!DTANTEpi$scluster %in% c("PT1", "PT2", "PT3")]
keep_ind <- !colnames(DTTEpi) %in% del_barcode
# Find markers for tumor and proximal tubule cells of each individual sample.
seurat2 <- DTTEpi[, keep_ind]
Idents(seurat2) <- "orig.ident"
if(!file.exists(individual_marker_file <- "epi_revise_individual_markers.rds")){
	individual_markers <- FindAllMarkers(seurat2, assay = "RNA", slot = "data", test.use = "MAST", 
		only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, return.thresh = 0.05)
	saveRDS(individual_markers, file = individual_marker_file)
}
saveRDS(seurat2, file = "epi_revise_seurat2_raw.rds")

# Make a CellDataSet and do some preprocessing.
cds_counts <- GetAssayData(seurat2, assay = "RNA", slot = "counts")
if(!inherits(cds_counts, "sparseMatrix")) cds_counts <- as(cds_counts, "sparseMatrix")
cds_pd <- as(seurat2@meta.data, "AnnotatedDataFrame")
cds_fd <- as(data.frame(row.names = rownames(cds_counts), gene_short_name = rownames(cds_counts)), "AnnotatedDataFrame")
cds <- newCellDataSet(
	cds_counts, phenoData = cds_pd, featureData = cds_fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
LOX_vec <- seurat2@assays$RNA@data["LOX", colnames(cds)]
if(anyNA(LOX_vec)) stop("Something went wrong", rep("!\n", 10))
cds$LOX <- LOX_vec
saveRDS(cds, file = "epi_revise_cds_raw.rds")

cds <- readRDS(file = "epi_revise_cds_raw.rds")
# Define up-regulated genes in ANT/TT1-TT/TT2-TT/TT1-PT/TT2-PT groups.
ident2group <- sapply(X = split(cds$TTCluster, f = cds$orig.ident), FUN = unique)
individual_markers <- readRDS(file = "epi_revise_individual_markers.rds")
levels(individual_markers$cluster) <- ident2group[levels(individual_markers$cluster)]
individual_markers_list <- split(individual_markers, f = sign(individual_markers$avg_log2FC))
group_markers_up_list <- lapply(X = split(individual_markers_list[["1"]], f = individual_markers_list[["1"]]$cluster), FUN = "[[", i = "gene")
# Ensure a selected gene for a group is expressed by >1 samples within that group.
group_markers_up_list <- lapply(X = group_markers_up_list, FUN = function(x){y <- table(x); names(y)[y > 1]})
# Eliminate genes shared between ANT, TT1 and TT2. Please refer to the methods section of our manuscript.
group_markers_up4 <- with(
	data = group_markers_up_list, 
	expr = unique(c(
		setdiff(ANT, c(`TT1-TT`, `TT1-PT`, `TT2-TT`, `TT2-PT`)), 
		setdiff(c(`TT1-TT`, `TT1-PT`), c(ANT, `TT2-TT`, `TT2-PT`)), 
		setdiff(c(`TT2-TT`, `TT2-PT`), c(ANT, `TT1-TT`, `TT1-PT`))
	))
)
saveRDS(group_markers_up4, file = "epi_revise_ordering_genes.rds")

save_diag_plot <- function(cds, suffix = "1st"){
	set.seed(123)
	p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = TRUE, cell_size = 0.5)
	p1 <- p1 + theme(legend.position = "bottom") + scale_colour_viridis_c(option = "inferno", name = "Pseudotime")
	p2 <- plot_cell_trajectory(cds, color_by = "TTCluster", show_branch_points = TRUE, cell_size = 0.5)
	p3 <- plot_cell_trajectory(cds, color_by = "State", show_branch_points = TRUE, cell_size = 0.5)
	p4 <- plot_cell_trajectory(cds, color_by = "orig.ident", show_branch_points = TRUE, cell_size = 0.5)
	p5 <- plot_cell_trajectory(cds, color_by = "LOX", show_branch_points = TRUE, cell_size = 0.5)
	p5 <- p5 + viridis::scale_color_viridis(option = "F", direction = -1)
	p5$layers[[2]]$aes_params$shape <- 16
	p5$layers[[2]]$aes_params$alpha <- 0.5
	p1r <- ggrastr::rasterise(p1, layers = c("Point"), dpi = 360, dev = "cairo", scale = 1)
	p <- patchwork::wrap_plots(p1r, p2, p3, p4, p5, nrow = 1)
	ggsave(p, filename = paste0("epi_revise_ordering_", suffix, ".png"), width = 32, height = 9)
	ggsave(p, filename = paste0("epi_revise_ordering_", suffix, ".pdf"), width = 32, height = 9)
}
# First attempt to order cells.
set.seed(123)
group_markers_up4 <- readRDS(file = "epi_revise_ordering_genes.rds")
cds <- setOrderingFilter(cds, unique(c(group_markers_up4)))
#cds <- reduceDimension(cds, norm_method = "vstExprs", verbose = TRUE)
cds <- reduceDimension(cds, verbose = TRUE)
cds <- orderCells(cds)
save_diag_plot(cds, suffix = "1st")
saveRDS(cds, file = "epi_revise_cds_1st.rds")

# Second attempt to order cells.
theoreticalNaiveState <- function(cds, ident, naiveIdent){
	if(length(unique(pData(cds)$State))){
		T0_counts <- table(cds$State, cds[[ident]])[, naiveIdent]
		return(as.numeric(names(T0_counts)[which.max(T0_counts)]))
	}else{
		return(1)
	}
}
currentNaiveState <- function(cds){
	x <- sapply(X = split(cds$Pseudotime, f = cds$State), FUN = mean, na.rm = TRUE)
	as.numeric(names(x)[which.min(x)])
}
theoretical_naive_state <- theoreticalNaiveState(cds, ident = "TTCluster", naiveIdent = "ANT")
if(currentNaiveState(cds) != theoretical_naive_state){
	set.seed(123)
	cds <- orderCells(cds, root_state = theoretical_naive_state)
	save_diag_plot(cds, suffix = "2nd")
	saveRDS(cds, file = "epi_revise_cds_2nd.rds")
}

# Smoothed expression curves and heatmaps along the pseudotime.
TT_color <- c(
	"TT1-TT" = "#ff7f0f", "TT1-PT" = "#ffcc9f", 
	"TT2-TT" = "#3cb7cc", "TT2-PT" = "#d5eff4", 
	"ANT" = "#87cfa4")
p <- ggplot(
	data = as.data.frame(pData(cds)[cds$TTCluster %in% c("ANT", "TT1-TT", "TT2-TT"), c("Pseudotime", "State", "TTCluster", "LOX")]), 
	aes(x = Pseudotime, y = LOX, colour = TTCluster, fill = TTCluster), 
	alpha = 0.05
) + geom_smooth(
	method = "loess"
)
p <- p + scale_color_manual(values = TT_color) + scale_fill_manual(values = TT_color)
p <- p + theme(aspect.ratio = 1)
ggsave(p, filename = "epi_revise_cds_curve.png", height = 5, width = 5)
ggsave(p, filename = "epi_revise_cds_curve.pdf", height = 5, width = 5)
# https://github.com/cole-trapnell-lab/monocle-release/issues/133
# need to plot at least two genes for plot_genes_branched_pseudotime to avoid the above error
filename <- "epi_revise_cds_heat"
gene_vec <- c("LOX", "CEBPB", "TGFB1", "PLOD2")
pdf(paste0(filename, ".pdf"), width = 8, height = 0.5 * length(gene_vec))
plot_genes_branched_heatmap(cds[gene_vec, ], num_clusters = 1, branch_point = 2, show_rownames = TRUE)
dev.off()
png::writePNG(pdftools::pdf_render_page(paste0(filename, ".pdf"), page = 1, dpi = 150), paste0(filename, ".png"))


## Cross study comparisons.
# Load the reference dataset.
seurat2 <- readRDS(file = "epi_revise_seurat2_raw.rds")
DefaultAssay(seurat2) <- "RNA"
cds <- readRDS(file = "epi_revise_cds_1st.rds")
# Load an integrated query dataset.
fuxkRDA <- function(file){args <- c(ls(), "args"); load(file); if(length(obj <- setdiff(ls(), args)) != 1L) stop(); get(obj);}
f <- "some.integrated.rda"
extern_data <- fuxkRDA(f)
DefaultAssay(extern_data) <- "RNA"

extern_data$dataset <- NA_character_
extern_data$dataset[grep("PD", extern_data$orig.ident)] <- "PMID36423636"
extern_data$dataset[which(extern_data$orig.ident %in% c("N1","N2","T1","T2","T3","T4","T5","T6","T7"))] <- "GSE210038"
extern_data$dataset[which(extern_data$orig.ident %in% c("GSM4735364_RCC1t", "GSM4735365_RCC1n"
,"GSM4735366_RCC2t", "GSM4735367_RCC2n", "GSM4735368_RCC3t", "GSM4735369_RCC3n"
,"GSM4735370_RCC4t", "GSM4735371_RCC4n", "GSM4735372_RCC5t", "GSM4735373_RCC5n"
, "GSM4735374_RCC6t", "GSM4735375_RCC7t"))] <- "GSE156632"
extern_data$dataset[which(extern_data$orig.ident %in% c("GSM7028034_RCC1", "GSM7028035_RCC2" 
,"GSM7028036_RCC3","GSM7028037_RCC4","GSM7028038_RCC4f","GSM7028039_RCC5","GSM7028040_RCC5t"))] <- "GSE224630"

# 
DDRT_mat <- t(reducedDimS(cds))
colnames(DDRT_mat) <- paste0("DDR_", 1:ncol(DDRT_mat))
if(!all(rownames(DDRT_mat) == colnames(seurat2))) stop("Something went wrong", rep("!\n", 10))
if(FALSE){
DDRT <- new("DimReduc", cell.embeddings = DDRT_mat, assay.used = "RNA", global = FALSE, key = "DDR_")
seurat2@reductions$DDRTree <- DDRT
seurat2@meta.data$Pseudotime <- cds$Pseudotime
seurat2@meta.data$State <- cds$State
p1 <- DimPlot(seurat2, reduction = "DDRTree", group.by = "State", shuffle = TRUE, raster = FALSE, pt.size = 0.05)
p2 <- FeaturePlot(seurat2, reduction = "DDRTree", features = c("LOX"), raster = FALSE, pt.size = 0.05)
ggsave(p1 + p2, filename = "epi_revise_seurat2_test.png", height = 6, width = 12)
saveRDS(seurat2, file = "epi_revise_seurat2_1st.rds")
}

# Seurat label transfer.
group_markers_up4 <- readRDS(file = "epi_revise_ordering_genes.rds")
TT.anchors <- FindTransferAnchors(
	reference = seurat2, query = extern_data, reduction = "rpca", features = group_markers_up4, 
	reference.assay = "RNA", query.assay = "RNA")
# Transfer DDRTree embedding.
predictions <- TransferData(anchorset = TT.anchors, refdata = t(DDRT_mat))
extern_DDRT_mat <- t(as.matrix(predictions@data))
colnames(extern_DDRT_mat) <- paste0("DDR_", 1:ncol(extern_DDRT_mat))
extern_DDRT <- new("DimReduc", cell.embeddings = extern_DDRT_mat, assay.used = "RNA", global = FALSE, key = "DDR_")
extern_data@reductions$DDRTree <- extern_DDRT
# Transfer Pseudotime.
predictions <- TransferData(anchorset = TT.anchors, refdata = t(data.matrix(pData(cds)[, "Pseudotime", drop = FALSE])))
extern_DDRT_time <- t(as.matrix(predictions@data))
extern_data@meta.data$Pseudotime <- extern_DDRT_time[, 1]

# 
p0 <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = TRUE, cell_size = 0.5)
p1 <- FeaturePlot(extern_data, reduction = "DDRTree", features = "Pseudotime", raster = TRUE, pt.size = 1.5)
p1 <- p1 + scale_colour_viridis_c(option = "inferno", name = "Pseudotime")
p2 <- DimPlot(extern_data, reduction = "DDRTree", group.by = "group", shuffle = TRUE, raster = TRUE, pt.size = 1.5)
p3 <- DimPlot(extern_data, reduction = "DDRTree", group.by = "dataset", shuffle = TRUE, raster = TRUE, pt.size = 1.5)
p4 <- FeaturePlot(extern_data, reduction = "DDRTree", features = c("LOX"), raster = TRUE, pt.size = 1.5)
p4 <- p4 + viridis::scale_color_viridis(option = "F", direction = -1)
p <- patchwork::wrap_plots(p1, p2, p3, p4, nrow = 1) & theme(aspect.ratio = 1, legend.position = "top")
for(i in 1:4) p[[i]]$layers <- c(p[[i]]$layers[[1]], p0$layers[[1]])
ggsave(p, filename = "epi_revise_extern_test.png", width = 17, height = 5)
ggsave(p, filename = "epi_revise_extern_test.pdf", width = 17, height = 5)
