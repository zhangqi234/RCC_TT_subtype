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

# Tested under R4.0. 
# May break easily as they are not well written and hijack too much...

# Draw a border around selected spots.
#' @param p is a ggplot object from \code{Seurat::SpatialDimPlot(object, ...)}
#' @param idents is the identity column to use in \code{object}
#' @param highlights is a vector of niches to be encircled within \code{idents}
#' @param color is the color of dots consisting the spatial border
addSpatialBorders <- function(
	p, 
	idents, 
	highlights, 
	color = "black", 
	...
){
	# Constants
	probs <- c(0, 1)
	max_neighbors <- 6L
	maxgap <- 1L
	
	# 'p' must inherit from 'ggplot'
	if(!is(p, "ggplot"))
		stop("addSpatialBorders: p has to be an object of ggplot")
	# Get the number of subplots if 'p' is a patchwork object
	p_idx <- 1L
	if(is(p, "patchwork")){
		p_idx <- seq_len(length(p$patches$plots) + 1L)
	}else{
		# Place 'p' in a list for convenience if it's a ggplot
		# Because we will access each subplot by 'p[[p_idx[i]]]'
		p <- list(p)
	}
	# Get the spatial layer in each subplot
	s_idx <- sapply(
		simplify = TRUE, 
		X = p_idx, 
		FUN = function(i){
			idx <- sapply(
				simplify = TRUE, 
				X = p[[i]]$layers, 
				FUN = function(x) is(x$geom, "GeomSpatial")
			)
			if(!any(idx)) stop(
				"addSpatialBorders: could not find the GeomSpatial layer ", 
				"in subplot ", i
			)
			idx_which <- which(idx)
			if(length(idx_which) > 1L) warning(
				"addSpatialBorders: 2 or more GeomSpatial layers ", 
				"were found in subplot ", i, ", ", 
				"but only the 1st one will be used"
			)
			return(idx_which)
		}
	)
	# Get the spot diameter on each slide
	diameter <- mapply(
		SIMPLIFY = TRUE, 
		p_idx, s_idx, 
		FUN = function(x_p, x_s) with(
			p[[x_p]]$layers[[x_s]]$geom_params$image@scale.factors, 
			lowres * fiducial
		)
	)
	diameter <- ceiling(maxgap * diameter)

	# Check validity of 'idents'
	if(missing(idents)){
		# If 'idents' is not supplied, data in each subplot should have an 'ident' column
		# Issue an error if this column could not be found
		if(!all(sapply(p_idx, function(i) "ident" %in% colnames(p[[i]]$data))))
			stop("addSpatialBorders: idents not provided")
	}else{
		# When 'idents' is supplied we overwrite the 'ident' column in each subplot
		override <- FALSE
		if(!is.factor(idents)) idents <- factor(idents)
		for(i in p_idx){
			if("ident" %in% colnames(p[[i]]$data)) override <- TRUE
			if(
				is.null(names(idents)) || 
				!all(rownames(p[[i]]$data) %in% names(idents))
			){
				stop("addSpatialBorders: incompatible idents")
			}
			p[[i]]$data$ident <- unname(idents[rownames(p[[i]]$data)])
		}
		# Issue a warning if any subplot data gets overwritten
		if(override) warning("addSpatialBorders: overwriting idents in the ggplot object")
	}
	
	# Encircle each ident respectively if 'highlights' not specified
	if(missing(highlights)) highlights <- as.list(levels(p$data$ident))
	if(!is.list(highlights))
		stop("addSpatialBorders: highlights has to be a list of clusters")

	# Create the spot matrix
	for(x_p in p_idx){
		# Get valid spot barcodes on this slide
		x_s <- s_idx[x_p]
		x_D <- diameter[x_p]
		barcodes <- intersect(
			rownames(p[[x_p]]$data), 
			rownames(p[[x_p]]$layers[[x_s]]$geom_params$image@coordinates)
		)
		# Generate meta.data for those spots
		meta_data <- cbind.data.frame(
			p[[x_p]]$data[barcodes, , drop = FALSE], 
			p[[x_p]]$layers[[x_s]]$geom_params$image@coordinates[barcodes, c("row", "col"), drop = FALSE]
		)
		levels(meta_data$ident) <- c(levels(meta_data$ident), "VIRTUAL")
		# Fix coordinates
		if(min(meta_data$row) == 0L) meta_data$row <- meta_data$row + 1L
		if(min(meta_data$col) == 0L) meta_data$col <- meta_data$col + 1L
		meta_data$imagerow <- max(meta_data$imagerow) + min(meta_data$imagerow) - meta_data$imagerow

		# Determine the 'ident' of neighbors surrounding each spot
		neighbors <- matrix(
			NA_integer_, 
			nrow = nrow(meta_data), 
			ncol = length(levels(meta_data$ident)), 
			dimnames = list(rownames(meta_data), levels(meta_data$ident))
		)
		for(i in rownames(meta_data)){
			a <- meta_data$imagerow - meta_data[i, "imagerow"]
			b <- meta_data$imagecol - meta_data[i, "imagecol"]
			d <- sqrt(a ^ 2 + b ^ 2)
			neighbors[i, ] <- as.vector(table(meta_data$ident[which(d <= x_D / maxgap)]))
			j <- meta_data[i, "ident", drop = TRUE]
			# Don't count the spot itself
			neighbors[i, j] <- neighbors[i, j] - 1L
		}
		# A spot has virtual neighbors if it's not surrounded by 6 other spots
		if(max(neighbors) < max_neighbors)
			warning("addSpatialBorders: check your data because all spots have <6 neighbors")
		if(max(neighbors) > max_neighbors)
			stop("addSpatialBorders: spots have >6 neighbors")
		neighbors[, "VIRTUAL"] <- neighbors[, "VIRTUAL"] + max_neighbors - rowSums(neighbors)
		# Count neighbors of different 'ident' for each spot
		n_self <- neighbors[(as.integer(meta_data$ident) - 1L) * nrow(neighbors) + seq_len(nrow(neighbors))]
		meta_data$n_neighbors <- max_neighbors - n_self

		# Determine the dimension of this slide
		row_limits <- quantile(meta_data$row, probs = probs)
		col_limits <- quantile(meta_data$col, probs = probs)
		# Create the spot matrix storing ident values
		slide <- matrix(nlevels(meta_data$ident), nrow = row_limits[2], ncol = col_limits[2])
		for(j in seq_len(nrow(meta_data))){
			slide[meta_data[j, "row"], meta_data[j, "col"]] <- meta_data[j, "ident"]
		}
		attributes(slide) <- c(attributes(slide), attributes(meta_data$ident))
		p[[x_p]]$spots.matrix <- slide
		p[[x_p]]$spots.neighbors <- neighbors
		p[[x_p]]$meta.data <- meta_data
	}

	# Add the borders to the subplots
	for(x_p in p_idx){
		x_s <- s_idx[x_p]
		x_D <- diameter[x_p]
		inds <- paste(p[[x_p]]$meta.data$row, p[[x_p]]$meta.data$col, sep = "\r")
		for(highlight in highlights){
			x <- integer(length(inds) * 2L)
			x[] <- NA_integer_
			y <- 0L
			border_spots <- data.frame(row = x, col = x)

			# Spots other than 'highlight' get NA
			spots <- p[[x_p]]$spots.matrix
			levels(spots)[which(!levels(spots) %in% highlight)] <- NA
			spots <- !is.na(spots)
			if(!any(spots)) next

			# Scan TRUE/FALSE borders by row (actually, cols in the image)
			for(i in seq_len(nrow(spots))){
				# A spot is in it's cluster if spots on both sides have the same ident
				# Here we define 'shade' as the identity of spots on both sides
				if(i == 1L){
					shade <- spots[2L, , drop = TRUE]
				}else if(i == nrow(spots)){
					shade <- spots[nrow(spots) - 1L, , drop = TRUE]
				}else{
					shade <- spots[i - 1L, , drop = TRUE] & spots[i + 1L, , drop = TRUE]
				}
				i_rle <- rle(spots[i, , drop = TRUE] | shade)

				# Skip when there is no 'TRUE' segment
				idx <- which(i_rle$values)
				if(length(idx) == 0L) next
				# Obtain the starting and ending coordinates of 'TRUE' segments
				i_end <- cumsum(i_rle$lengths)
				i_start <- c(1L, i_end[-length(i_end)] + 1L)
				i_end <- i_end[idx]
				i_start <- i_start[idx]

				# If a segment starts or ends due to a shaded spot, shrink it
				idx <- which(shade & ! spots[i, , drop = TRUE])
				i_rm <- which(i_end %in% idx)
				i_end[i_rm] <- i_end[i_rm] - 1L
				i_rm <- which(i_start %in% idx)
				i_start[i_rm] <- i_start[i_rm] + 1L

				# Ensure ending pos is greater than starting pos
				idx <- which(i_start <= i_end)
				if(length(idx) == 0L) next
				i_end <- i_end[idx]
				i_start <- i_start[idx]

				# Write the starting/ending pos of all 'TRUE' segments
				y <- seq(from = y + 1L, to = y + length(idx) * 2L, by = 1L)
				border_spots[y, "row"] <- i
				border_spots[y, "col"] <- c(i_start, i_end)
				y <- y[length(y)]
			}
			
			# Remove duplicated border spots
			border_spots <- border_spots[which(!is.na(border_spots$row)), ]
			r_name <- do.call(paste, c(border_spots, sep = "\r"))
			idx <- which(!duplicated(r_name))
			border_spots <- border_spots[idx, ]
			rownames(border_spots) <- r_name[idx]
			border_spots <- border_spots[which(rownames(border_spots) %in% inds), ]
			border_spots <- p[[x_p]]$meta.data[match(rownames(border_spots), inds), ]

			# Add border points to plot
			p[[x_p]] <- p[[x_p]] + geom_point(
				data = border_spots[, c("imagecol", "imagerow")], 
				mapping = aes(x = imagecol, y = imagerow), 
				inherit.aes = FALSE, 
				...
			)
		}
	}
	if(is(p, "patchwork")) return(p)
	p[[1L]]
}

getSpatialBorders <- function(p){
	# 'p' must inherit from 'ggplot'
	if(!is(p, "ggplot"))
		stop("getSpatialBorders: p has to be an object of ggplot")
	# Get the number of subplots if 'p' is a patchwork object
	p_idx <- 1L
	if(is(p, "patchwork")){
		p_idx <- seq_len(length(p$patches$plots) + 1L)
	}else{
		# Place 'p' in a list for convenience if it's a ggplot
		# Because we will access each subplot by 'p[[p_idx[i]]]'
		p <- list(p)
	}

	# Get the border point layer in each subplot
	s_idx <- sapply(
		simplify = TRUE, 
		X = p_idx, 
		FUN = function(i){
			idx <- sapply(
				simplify = TRUE, 
				X = p[[i]]$layers, 
				FUN = function(x) is(x$geom, "GeomPoint")
			)
			if(!any(idx)) stop(
				"getSpatialBorders: could not find the GeomPoint layer ", 
				"in subplot ", i
			)
			idx_which <- which(idx)
			if(length(idx_which) > 1L) warning(
				"getSpatialBorders: 2 or more GeomPoint layers ", 
				"were found in subplot ", i, ", ", 
				"but only the 1st one will be used"
			)
			return(idx_which)
		}
	)

	sapply(
		simplify = TRUE, 
		X = p_idx, 
		FUN = function(x_p){
			x_s <- s_idx[x_p]
			x_i <- p[[x_p]]$layers[[x_s - 1L]]$geom_params$image@key
			barcodes <- rownames(p[[x_p]]$layers[[x_s]]$data)
			v <- list(barcodes)
			names(v) <- x_i
			return(v)
		}
	)
}

getSpatialNeighbors <- function(p){
	# 'p' must inherit from 'ggplot'
	if(!is(p, "ggplot"))
		stop("getSpatialNeighbors: p has to be an object of ggplot")
	# Get the number of subplots if 'p' is a patchwork object
	p_idx <- 1L
	if(is(p, "patchwork")){
		p_idx <- seq_len(length(p$patches$plots) + 1L)
	}else{
		# Place 'p' in a list for convenience if it's a ggplot
		# Because we will access each subplot by 'p[[p_idx[i]]]'
		p <- list(p)
	}

	# Get the spatial layer in each subplot
	s_idx <- sapply(
		simplify = TRUE, 
		X = p_idx, 
		FUN = function(i){
			idx <- sapply(
				simplify = TRUE, 
				X = p[[i]]$layers, 
				FUN = function(x) is(x$geom, "GeomSpatial")
			)
			if(!any(idx)) stop(
				"getSpatialNeighbors: could not find the GeomSpatial layer ", 
				"in subplot ", i
			)
			idx_which <- which(idx)
			if(length(idx_which) > 1L) warning(
				"getSpatialNeighbors: 2 or more GeomSpatial layers ", 
				"were found in subplot ", i, ", ", 
				"but only the 1st one will be used"
			)
			return(idx_which)
		}
	)

	sapply(
		simplify = TRUE, 
		X = p_idx, 
		FUN = function(x_p){
			x_s <- s_idx[x_p]
			x_i <- p[[x_p]]$layers[[x_s]]$geom_params$image@key
			df <- p[[x_p]]$spots.neighbors
			df <- df[, -which(colnames(df) == "VIRTUAL")]
			df <- list(df)
			names(df) <- x_i
			return(df)
		}
	)
}

# A improved version of \code{Seurat::SpatialDimPlot} which allows more 
# controls and customizations.
SpatialDimPlot <- function(
	object, 
	..., 
	expand = expansion(mult = c(0.02, 0.02), add = c(0, 0)), 
	greyscale = FALSE, 
	pt.size = 2, 
	pt.shape = 16, 
	rasterize = FALSE
){
	p <- Seurat::SpatialDimPlot(object, ...)
	for(i in seq.int(length(p))){
		# Backup layers (so that they won't be influenced by the news cale)
		l_st <- p[[i]]$layers[[1L]]
		l_list <- list()
		if(length(p[[i]]$layers) > 1) l_list <- p[[i]]$layers[-1]
		p[[i]]$layers <- list()
		# Get image object
		idx <- sapply(
			X = object@images, 
			FUN = function(img) img@key == l_st$geom_params$image@key
		)
		img <- object@images[[which(idx)]]
		for(j in ls(p[[i]]$plot_env)){
			get(j, envir = p[[i]]$plot_env)
			rm(list = j, envir = p[[i]]$plot_env)
		}
		# Get center (flip y)
		if(is(img, "VisiumV1")){
			pt_df <- img@coordinates
		}else if(is(img, "VisiumV2")){
			pt_df <- data.frame(
				imagerow = img@boundaries[[1]]@coords[, "x"], 
				imagecol = img@boundaries[[1]]@coords[, "y"]
			)
		}else{
			stop("unknown img")
		}
		pt_df$imagerow <- pt_df$imagerow * img@scale.factors$lowres
		pt_df$imagecol <- pt_df$imagecol * img@scale.factors$lowres
		y_cen <- sum(range(pt_df$imagerow))
		if(is(img, "VisiumV1")){
			p[[i]]$data$imagerow <- y_cen - p[[i]]$data$imagerow
		}else{
			p[[i]]$data$x <- y_cen - p[[i]]$data$x
		}

		# Plot the image
		col_mat <- img@image
		dimnames(col_mat) <- list(NULL, NULL, c('r', 'g', 'b'))
		if(greyscale) for(a in seq.int(nrow(col_mat))) for(b in seq.int(ncol(col_mat))){
			col_mat[a, b, ] <- sum(col_mat[a, b, ] * c(0.2989, 0.5870, 0.1140))
		}
		col_df <- reshape2::melt(col_mat[, , 1], value.name = "r")
		col_df$g <- as.vector(col_mat[, , 2])
		col_df$b <- as.vector(col_mat[, , 3])
		col_df$Var1 <- y_cen - col_df$Var1
		p[[i]] <- p[[i]] + geom_raster(
			data = col_df, 
			mapping = aes(
				x = Var2, 
				y = Var1, 
				fill = rgb(r, g, b, alpha = 1, maxColorValue = 1)
			), 
			alpha = l_st$geom_params$image.alpha, 
			hjust = 0, vjust = 0, 
			inherit.aes = FALSE
		)
		# Backup and restore the 'fill' scale
		idx <- sapply(X = p[[i]]$scales$scales, function(s) s$aesthetics == "fill")
		s_backup <- NULL
		if(any(idx))
			s_backup <- p[[i]]$scales$scales[[which(idx)[1]]]
		p[[i]] <- p[[i]] + scale_fill_identity()
		p[[i]] <- p[[i]] + ggnewscale::new_scale_fill()
		if(is(s_backup, "Scale"))
			p[[i]]$scales$scales[[length(p[[i]]$scales$scales) + 1L]] <- s_backup
		p[[i]]$mapping$fill <- p[[i]]$mapping$fill_new
		p[[i]]$mapping$fill_new <- NULL

		# Add spots
		if(is.character(pt.shape)){
			p[[i]] <- p[[i]] + ggstar::geom_star(
				starshape = pt.shape, 
				starstroke = 0, 
				angle = 90, 
				size = pt.size, 
				alpha = l_st$aes_params$alpha, 
				colour = NA
			)
		}else{
			p[[i]] <- p[[i]] + geom_point(
				shape = pt.shape, 
				size = pt.size, 
				alpha = l_st$aes_params$alpha, 
				fill = NA
			)
		}
		if(!is.null(p[[i]]$mapping$alpha))
			# Let it inherit the default alpha mapping
			p[[i]]$layers[[length(p[[i]]$layers)]]$aes_params$alpha <- NULL
		# Copy the 'colour' from 'fill' scale
		p[[i]]$mapping$colour <- p[[i]]$mapping$fill
		if(is(s_backup, "Scale")){
			s_colour <- structure(new.env(), class = class(s_backup))
			for(o in ls(s_backup)) assign(o, get(o, envir = s_backup), envir = s_colour)
			s_colour$aesthetics <- "colour"
			p[[i]]$scales$scales[[length(p[[i]]$scales$scales) + 1L]] <- s_colour
		}

		# Add backuped layers
		p[[i]]$layers <- c(p[[i]]$layers, l_list)

		# Restrict axis
		x_min <- floor(min(pt_df$imagecol))
		x_max <- ceiling(max(pt_df$imagecol))
		y_min <- floor(min(pt_df$imagerow))
		y_max <- ceiling(max(pt_df$imagerow))
		x_add_l <- (x_max - x_min) * expand[1] + expand[2]
		x_add_r <- (x_max - x_min) * expand[3] + expand[4]
		y_add_b <- (x_max - x_min) * expand[1] + expand[2]
		y_add_t <- (x_max - x_min) * expand[3] + expand[4]
		p[[i]] <- p[[i]] + scale_y_continuous(
			limits = c(
				y_cen - (y_max + y_add_t), 
				y_cen - (y_min - y_add_b)
			)
		)
		p[[i]] <- p[[i]] + scale_x_continuous(
			limits = c(
				x_min - x_add_l, 
				x_max + x_add_r
			)
		)
		p[[i]]$coordinates$expand <- FALSE
		p[[i]] <- p[[i]] + theme(
			panel.border = element_rect(colour = "black", fill = "transparent")
		)

		# Rasterize
		if(rasterize) p[[i]] <- ggrastr::rasterise(p[[i]], layers = c("Star", "Point"))
	}
	return(p)
}
