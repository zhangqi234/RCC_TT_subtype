'
MIT License

Copyright (c) 2024 Tan Yezhen, <11510390@mail.sustech.edu.cn>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'

# Dependencies
library("ggplot2")
stopifnot(all(c(
	"patchwork", 
	"viridis", 
	"ggforce"
) %in% installed.packages()))

# This function takes Visium spots and their identities, base on which it
# computes geographically smoothed identity scores
#' @param object is the Seurat object with Visium V1 data
#' @param idents is the meta.data column or a logical vector giving the 
# identity of individual spots (e.g. for a 'tumor-or-not' identity, an identity
# value of 1 specifies a tumor spot while zeros stand for none-tumor spots)
#' @param d is the half-width of the hexagonal sliding title
#' @param compare is a Boolean value indicating whether heatmaps of the tissue
# section should be drawn to compare the smoothed identity values with the 
# original ones
#' @note 
getGeoScores <- function(
	object, 
	idents, 
	d = Inf, 
	compare = FALSE
){
	# Sanity check
	if(length(object@images) == 0L)
		stop("spatialRegionScore: object does not contain spatial information")
	if(length(idents) == 1L)
		idents <- object@meta.data[, idents, drop = TRUE]
	if(length(idents) != ncol(object))
		stop("spatialRegionScore: idents does not match object")
	if(!is.logical(idents))
		idents <- as.logical(idents)

	# The return value is a named numeric vector
	geo_score <- numeric(length = ncol(object))
	names(geo_score) <- colnames(object)
	geo_score[] <- NA_real_

	# For each tissue section, smooth each spot with sliding tiles
	# (like a sliding window, but on the 2D surface)
	for(i in names(object@images)){
		# Get valid barcodes on this tissue section
		barcodes <- intersect(
			rownames(object@images[[i]]@coordinates), 
			names(geo_score)
		)
		mappings <- match(barcodes, names(geo_score))

		# Generate meta.data for those spots
		meta_data <- cbind.data.frame(
			object@images[[i]]@coordinates[barcodes, ], 
			ident = idents[mappings]
		)

		# Ensure coordinates start from (1, 1)
		if(min(meta_data$row) != 1L)
			meta_data$row <- meta_data$row - min(meta_data$row) + 1L
		if(any(meta_data$col < 1L))
			meta_data$col <- meta_data$col - min(meta_data$col) + 1L

		# Check the half-width of the sliding tile
		if(!is.finite(d)) d <- min(
			max(meta_data$row), 
			max(meta_data$col) %/% 2L
		) %/% 2L - 1L
		# Ensure the half-width >= 0L
		if(d < 0L) d <- 0L

		# Add padding to coordinates
		meta_data$row <- meta_data$row + d
		meta_data$col <- meta_data$col + d * 2L

		# Create a weight matrix for the sliding tile
		w_matrix <- getGeoWeightsT(d)

		# Reconstruct the spatial matrix
		s_matrix_s <- s_matrix <- s_name <- matrix(
			NA_integer_, 
			nrow = max(meta_data$row) + d, 
			ncol = max(meta_data$col) + d * 2L
		)
		mode(s_name) <- "character"
		for(s in seq.int(nrow(meta_data))){
			x <- meta_data[s, "row"]
			y <- meta_data[s, "col"]
			s_matrix[x, y] <- meta_data[s, "ident"]
			s_name[x, y]   <- rownames(meta_data)[s]
		}

		# Moving tile averaging of the spatial matrix with the weight matrix
		for(x in seq.int(nrow(s_matrix_s))){
			for(y in seq.int(ncol(s_matrix_s))){
				if(is.na(s <- s_name[x, y])) next
				tmp <- s_matrix[
					seq(from = x - d, to = x + d, by = 1L), 
					seq(from = y - d * 2L, to = y + d * 2L, by = 1L)
				]
				# Uncertain about the identities of the spots with NAs (
				# whether they are tumor spots or not?).
				# So the safest way is to exclude those spots prior to weighted
				# averaging.
				idx <- !is.na(tmp)
				geo_score[s] <- s_matrix_s[x, y] <- sum(
					w_matrix[idx] * tmp[idx] / sum(w_matrix[idx]), 
					na.rm = TRUE
				)
			}
		}

		# If requested, plot the original and smoothed geographic matrix
		if(compare){
			to_plot <- cbind(
				t(s_matrix[rev(seq.int(nrow(s_matrix))), ]), 
				t(s_matrix_s[rev(seq.int(nrow(s_matrix_s))), ])
			)
			heatmap(to_plot, Rowv = NA, Colv = NA, scale = "none", asp = 2/3)
		}
	}
	# 'geo_score' is the smoothed density of TRUEs in 'ident'
	# To obtain the distance to TRUEs:
	return(1 - geo_score)
}

# This function creates an empty geographic matrix
emptyGeoRegions <- function(n_col){
	s_matrix <- matrix(
		NA_integer_, 
		nrow = n_col * 2L - 1L, 
		ncol = n_col
	)
	s_vec1 <- logical(length = nrow(s_matrix))
	s_vec2 <- s_vec1
	s_vec1[seq(from = 1L, to = length(s_vec1), by = 2L)] <- NA
	s_vec2[seq(from = 2L, to = length(s_vec2), by = 2L)] <- NA
	if((n_col %/% 2L) %% 2L == 0L){
		s_matrix[, seq(from = 2L, to = ncol(s_matrix), by = 2L)] <- s_vec1
		s_matrix[, seq(from = 1L, to = ncol(s_matrix), by = 2L)] <- s_vec2
	}else{
		s_matrix[, seq(from = 2L, to = ncol(s_matrix), by = 2L)] <- s_vec2
		s_matrix[, seq(from = 1L, to = ncol(s_matrix), by = 2L)] <- s_vec1
	}
	return(s_matrix)
}

# This function demonstrates how a tissue section is segmented based on the 
# identities of spatial spots
#' @param d is the half-width of the hexagonal sliding title
#' @param size is size of the spots in the resulting plot
#' @note If 'p' is the resulting plot, levels and breaks for the identity 
# values could be obtained via p$levels and p$breaks
delimitGeoRegions <- function(
	d, 
	size = 1, 
	linewidth = 1
){
	# Weights of each layer
	#d <- d + 1L
	#i <- seq.int(d)
	#w <- (i - 1L) * 6L * (d - i + 1L)
	#w[1L] <- d
	# Weight matrix
	w_matrix <- getGeoWeights(d)
	w_matrix <- w_matrix / sum(w_matrix)

	# Tumor core
	w_binary <- emptyGeoRegions(ncol(w_matrix))
	w_binary[w_matrix != 0] <- 1L

	# If there are d layers, generate a minimal matrix with a tumor core
	s_tmp <- emptyGeoRegions(ncol(w_matrix) * 3L)
	s_row <- (nrow(w_binary) + 1L):(nrow(w_binary) * 2L) + 1L
	s_col <- (ncol(w_binary) + 1L):(ncol(w_binary) * 2L)
	s_tmp[s_row, s_col] <- w_binary

	# Wrap this matrix with NAs
	s_bg <- matrix(
		NA_integer_, 
		nrow = ncol(w_matrix) * 10L - 1L, 
		ncol = ncol(w_matrix) * 5L
	)
	bg_row <- (nrow(w_binary) + 1L):(nrow(w_binary) + nrow(s_tmp)) + 1L
	bg_col <- (ncol(w_binary) + 1L):(ncol(w_binary) + ncol(s_tmp))
	s_bg[bg_row, bg_col] <- s_tmp
	s_matrix_s <- s_matrix <- s_bg

	# Moving tile averaging
	for(y in seq.int(nrow(s_matrix_s))){
		for(x in seq.int(ncol(s_matrix_s))){
			if(is.na(s_matrix[y, x])) next
			tmp <- w_matrix * s_matrix[
				seq(from = y - d * 2L, to = y + d * 2L, by = 1L), 
				seq(from = x - d, to = x + d, by = 1L)
			]
			s_matrix_s[y, x] <- sum(tmp, na.rm = TRUE) 
		}
	}

	# 
	w_matrix[w_matrix == 0] <- NA
	s_matrix <- 1 - s_matrix[bg_row, bg_col]
	s_matrix_s <- 1 - s_matrix_s[bg_row, bg_col]
	bks_raw <- round(as.vector(s_matrix_s), digits = 3L)
	bks_raw <- na.omit(unique(bks_raw))
	bks_raw <- sort(bks_raw, decreasing = FALSE)
	# Frontier of tumor niche
	front_score <- max(s_matrix_s[s_row, s_col][w_matrix != 0], na.rm = TRUE)
	front_score <- round(front_score, digits = 3L)
	# Peri-tumor niche
	s_tmp <- s_matrix_s
	s_tmp[s_row, s_col][w_matrix != 0] <- Inf
	peri_score <- min(s_tmp, na.rm = TRUE)
	peri_score <- round(peri_score, digits = 3L)
	# Verify front and peri position scores
	peri_pos <- which(bks_raw == peri_score)
	front_pos <- which(bks_raw == front_score)
	if(peri_pos == front_pos){
		warning("Selecting a large 'd' leads to nearly identical peri- and front-niche scores")
	}else if((peri_pos - 1L) != front_pos){
		warning("Selecting a large 'd' is discouraged")
		peri_score <- front_score
	}
	# New breaks
	bks <- c(inner = 0, front = front_score, peri = peri_score, distal = 1)
	# Plot this transformation
	p_pre_trans <- p_post_trans <- p_weights <- ggplot() + geom_point(
		data = na.omit(reshape2::melt(s_matrix)), 
		mapping = aes(
			x = Var2, 
			y = Var1, 
			colour = value
		), 
		shape = 16, 
		size = size
	)
	layer_tmp <- serialize(object = p_post_trans$layers[[1]], connection = NULL)
	p_post_trans$layers[[1]] <- unserialize(connection = layer_tmp, refhook = NULL)
	p_post_trans$layers[[1]]$data <- na.omit(reshape2::melt(s_matrix_s))
	p_weights$layers[[1]] <- unserialize(connection = layer_tmp, refhook = NULL)
	p_weights$layers[[1]]$data <- na.omit(reshape2::melt(w_matrix))
	p <- patchwork::wrap_plots(p_pre_trans, p_post_trans, nrow = 1, guides = "collect")
	p <- p & viridis::scale_colour_viridis(
		option = "F", 
		direction = -1, 
		limits = c(-0.01, 1.01), 
		breaks = bks, 
		labels = names(bks), 
		name = "Relative distance to tumor niche"
	)
	p <- p & theme_void() & theme(aspect.ratio = 1, panel.border = element_rect())
	p$breaks <- bks
	p$levels <- bks_raw
	# Encircle the tumor nest
	p <- p & ggforce::geom_ellipse(
		mapping = aes(
			x0 = ceiling(ncol(s_matrix) / 2L), 
			y0 = ceiling(nrow(s_matrix) / 2L), 
			a = d + 0.5, 
			b = d * 2 + 1, 
			angle = 0
		), 
		colour = "red", 
		linewidth = linewidth, 
		linetype= "dashed"
	)
	return(p)
}

# This function creates the weight matrix for a hexagonal sliding title with 
# half-width d
getGeoWeights <- function(d){
	w_matrix <- matrix(0L, nrow = 4L * d + 1L, ncol = 2L * d + 1L)
	w_matrix[2L * d + 1L, d + 1L] <- d + 1L
	for(layer in seq.int(d)){
		y_min <- layer
		y_max <- 2L * d + 1L - y_min + 1L
		x_min <- layer * 2L - 1L
		x_max <- 4L * d + 1L - x_min + 1L
		y <- d + 1L
		x <- x_min
		while(y != y_max){
			w_matrix[x, y] <- layer
			x <- x + 1L
			y <- y + 1L
		}
		for(a in seq.int(d - layer + 1L)){
			w_matrix[x, y] <- layer
			x <- x + 2L
		}
		while(x != x_max){
			w_matrix[x, y] <- layer
			x <- x + 1L
			y <- y - 1L
		}
		while(y != y_min){
			w_matrix[x, y] <- layer
			x <- x - 1L
			y <- y - 1L
		}
		for(a in seq.int(d - layer + 1L)){
			w_matrix[x, y] <- layer
			x <- x - 2L
		}
		while(x != x_min){
			w_matrix[x, y] <- layer
			x <- x - 1L
			y <- y + 1L
		}
	}
	return(w_matrix)
}

# Create the weight matrix, but transposed
getGeoWeightsT <- function(d){
	w_matrix <- matrix(0L, ncol = 4L * d + 1L, nrow = 2L * d + 1L)
	w_matrix[d + 1L, 2L * d + 1L] <- d + 1L
	for(layer in seq.int(d)){
		x_min <- layer
		x_max <- 2L * d + 1L - x_min + 1L
		y_min <- layer * 2L - 1L
		y_max <- 4L * d + 1L - y_min + 1L
		x <- d + 1L
		y <- y_min
		while(x != x_max){
			w_matrix[x, y] <- layer
			x <- x + 1L
			y <- y + 1L
		}
		for(a in seq.int(d - layer + 1L)){
			w_matrix[x, y] <- layer
			y <- y + 2L
		}
		while(y != y_max){
			w_matrix[x, y] <- layer
			x <- x - 1L
			y <- y + 1L
		}
		while(x != x_min){
			w_matrix[x, y] <- layer
			x <- x - 1L
			y <- y - 1L
		}
		for(a in seq.int(d - layer + 1L)){
			w_matrix[x, y] <- layer
			y <- y - 2L
		}
		while(y != y_min){
			w_matrix[x, y] <- layer
			x <- x + 1L
			y <- y - 1L
		}
	}
	return(w_matrix)
}
