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
stopifnot(all(c(
	"collapse"
) %in% installed.packages()))
library("GenomicRanges")
library("SummarizedExperiment")

#' Get the genomic ranges of chromosome arms from cytoband coordinates
#' @param cytoband A GenomicRanges, a data.frame or the path to a cytoband file
# chr_arms <- cytobandToArm("ucsc.hg19.cytoband.txt")
cytobandToArm <- function(cytoband){
	# Read cytoband
	if(is.character(cytoband))
		cytoband <- read.delim(file = cytoband[1], header = TRUE)
	if(is.data.frame(cytoband)){
		colnames(cytoband) <- c("chr", "start", "end", "name", "gieStain")
		cytoband$start <- cytoband$start + 1L
		cytoband <- as(cytoband, "GRanges")
	}else if(!is(cytoband, "GenomicRanges")) stop(
		"cytobandToArm: input should be a GenomicRanges, ", 
		"a data.frame or the path to a cytoband file"
	)
	# Obtain the range of chromosome arms
	arm_by_chr <- lapply(
		X = split(x = cytoband, f = seqnames(cytoband)), 
		FUN = function(cc){
			# Partition (p, q, cen) of cytobands
			cyto_id <- gsub(
				pattern = "^([pq]).*", 
				replacement = "\\1", 
				x = cc$name
			)
			cyto_id[grep(pattern = "cen", x = cc$gieStain)] <- "cen"
			if(!all(cyto_id %in% c("p", "q", "cen")))
				stop("cytobandToArm: cytoband identification error")
			# Reduce to arms
			arms <- sapply(X = split(x = cc, f = cyto_id), range)
			arms <- unlist(as(arms, "GRangesList"), use.names = TRUE)
			# Assign names
			arms$type <- names(arms)
			names(arms) <- paste0(seqnames(arms), ":", arms$type)
			return(arms)
		}
	)
	sort(unlist(as(arm_by_chr, "GRangesList"), use.names = FALSE))
}

#' Get the baseline of some values (e.g. copy numbers)
#' @param object is an instance of class SummarizedExperiment
#' @param assay is where the baseline value will be estimated
#' @param method A numeric value, or one of "mode", "mean", "median", 
#' "weightedMode", "weightedMean", and "weightedMedian". It defines the values 
#' in \code{assay} at neutral copy numbers. Default is \code{"mode"}
getBaseline <- function(
	object, 
	assay = "segmented", 
	method = "mode"
){
	# Fetch assay data and bin annotations
	tmp <- assay(object, assay)
	bin <- rowRanges(object)
	# Check methods
	if(is.character(method)){
		method <- match.arg(
			arg = method, 
			choices = c(
				"mode", "mean", "median", 
				"weightedMode", "weightedMedian", "weightedMean"
			), 
			several.ok = FALSE
		)
		return(switch(
			EXPR = method, 
			mode = 
				collapse::fmode(x = tmp, na.rm = TRUE, w = NULL), 
			weightedMode = 
				collapse::fmode(x = tmp, na.rm = TRUE, w = width(bin)), 
			mean = 
				collapse::fmean(x = tmp, na.rm = TRUE, w = NULL), 
			weightedMean = 
				collapse::fmean(x = tmp, na.rm = TRUE, w = width(bin)), 
			median = {
				mode(tmp) <- "numeric"
				collapse::fmedian(x = tmp, na.rm = TRUE, w = NULL)
			}, 
			weightedMedian = {
				mode(tmp) <- "numeric"
				collapse::fmedian(x = tmp, na.rm = TRUE, w = width(bin))
			}, 
			stop("getBaseline: unknown method ", method)
		))
	}else if(is.numeric(method)){
		return(method)
	}else{
		stop("getBaseline: unsupported 'method'")
	}
}

#' Chromosomal instability assessment
#' @param object A \code{QDNAseqSignals} or \code{RangedSummarizedExperiment}
#' object
#' @param assay Assay to use. Default is \code{"segmented"}
#' @param baseline Refer to \code{method} in \code{getBaseline}
#' @param ta Threshold for copy number amplifications. Values in \code{assay}
#' above \code{baseline} + this positive value are considered amplified. 
#' Default is \code{0.1}.
#' @param td Threshold for copy number deletions. Values in \code{assay} below
#' \code{baseline} - this positive value are considered deleted. Default is 
#' \code{0.1}.
#' @note 
genomeInstability <- function(
	object, 
	assay = "segmented", 
	baseline = "mode", 
	ta = 0.1, 
	td = 0.1
){
	# Fetch assay data and bin annotations
	tmp <- assay(object, assay)
	bin <- rowRanges(object)
	if(is.null(tmp))
		stop("genomeInstability: assay not found.")
	if(is.null(bin))
		stop("genomeInstability: genomic ranges not found.")
	stopifnot(isDisjoint(bin))
	# Determine bins to use
	if(is.logical(bin$use)){
		tmp <- tmp[bin$use, , drop = FALSE]
		bin <- bin[bin$use]
	}else{
		message(
			"genomeInstability: 'use' in the bin annotation ", 
			"is not logical. Using all."
		)
	}
	if(any(is.na(tmp)))
		message("genomeInstability: the assay contains missing values.")
	non_autosomes <- c(
		"23", "24", "25", 
		"X", "Y", "M", 
		"chr23", "chr24", "chr25", 
		"chrX", "chrY", "chrM"
	)
	if(any(non_autosomes %in% seqlevelsInUse(bin)))
		warning("genomeInstability: are you using sex or mt chromosomes?")

	# Get the baseline of each column
	blines <- getBaseline(object = object, assay = assay, method = baseline)
	# Clean
	rm(object); invisible(gc())
	# Subtract the baseline to get the residual
	tmp <- sweep(tmp, MARGIN = 2, STATS = blines, FUN = "-")
	# Determine abnormal bins
	tmp <- (tmp > ta) | (tmp < -td)

	# Now split bins according to chromosomes
	sen <- as.factor(seqnames(bin))
	tmp <- lapply(
		X = split(x = tmp, f = sen, drop = TRUE), 
		FUN = matrix, 
		ncol = ncol(tmp), 
		byrow = FALSE
	)
	bin <- split(x = bin, f = sen, drop = TRUE)

	# Fraction instable region at each chromosome
	ins_per_chr <- mapply(
		SIMPLIFY = FALSE, 
		tmp, bin, 
		FUN = function(value, region){
			len <- width(region)
			len_frac <- len / sum(len)
			return(colSums(value * len_frac, na.rm = TRUE))
		}
	)
	# Clean
	rm(tmp); invisible(gc())
	rm(bin); invisible(gc())
	do.call(cbind, ins_per_chr)
}

# Instability assessment at given genomic regions
#' @param object Refer to genomeInstability
#' @param region A \code{GenomicRanges} where fraction of amplifications and
#' deletions will be assessed
#' @param assay Refer to genomeInstability
#' @param baseline Refer to \code{method} in \code{getBaseline}
#' @param ta Refer to genomeInstability
#' @param td Refer to genomeInstability
#' @param minoverlap Minimal fraction of each region covered by \code{object}. 
#' Regions with coverage less than this value get \code{NA}
#' @return A list with two matrices, pct.amp and pct.del, which are the 
#' amplified and deleted fractions of \code{region} covered by \code{object}. 
#' The per-region coverage is stored in the \code{pct.overlap} attribute.
#' @note 
regionAneuploidy <- function(
	object, 
	region, 
	assay = "segmented", 
	baseline = "mode", 
	ta = 0.1, 
	td = 0.1, 
	minoverlap = 0.9
){
	# Fetch assay data and bin annotations
	tmp <- assay(object, assay)
	bin <- rowRanges(object)
	if(is.null(tmp))
		stop("regionAneuploidy: assay not found.")
	if(is.null(bin))
		stop("regionAneuploidy: genomic ranges not found.")
	stopifnot(isDisjoint(bin))
	# Determine bin to use
	if(is.logical(bin$use)){
		tmp <- tmp[bin$use, , drop = FALSE]
		bin <- bin[bin$use]
	}else{
		message(
			"regionAneuploidy: 'use' in the bin annotation ", 
			"is not logical. Using all."
		)
	}
	if(any(is.na(tmp)))
		message("regionAneuploidy: the assay contains missing values.")

	# Get the baseline of each column
	blines <- getBaseline(object = object, assay = assay, method = baseline)
	# Clean
	rm(object); invisible(gc())
	# Subtract the baseline to get the residual
	tmp <- sweep(tmp, MARGIN = 2, STATS = blines, FUN = "-")

	# 'pct.overlap' is the fraction of 'region' covered by 'bin'
	pct.overlap <- numeric()
	length(pct.overlap) <- length(region)
	pct.del <- pct.amp <- matrix(
		NA_real_, 
		nrow = length(region), 
		ncol = ncol(tmp), 
		dimnames = list(names(region), colnames(tmp))
	)
	i <- 1
	for(i in seq_along(region)){
		r <- region[i]
		o <- findOverlaps(query = r, subject = bin)
		w <- width(pintersect(r[queryHits(o)], bin[subjectHits(o)]))
		w_len <- sum(w)
		# Percentage of 'r' that overlaps 'bins'
		f <- w_len / width(r)
		# Assign this to 'pct.overlap'
		pct.overlap[i] <- f
		# Return NAs if pct.overlap is below some threshold 'minoverlap'
		if(f < minoverlap) next
		# The instability index is the altered fraction of the intersection
		tmp.resid <- tmp[subjectHits(o), , drop = FALSE]
		pct.del[i, ] <- colSums(w * (tmp.resid < -td), na.rm = TRUE) / w_len
		pct.amp[i, ] <- colSums(w * (tmp.resid >  ta), na.rm = TRUE) / w_len
	}
	intersected <- list(pct.del = pct.del, pct.amp = pct.amp)
	# Make the 'pct.overlap' attribute
	attr(intersected, "pct.overlap") <- pct.overlap

	# Clean
	rm(tmp); invisible(gc())
	rm(bin); invisible(gc())
	return(intersected)
}

# Call gained or lost regions
#' @param x is the result from regionAneuploidy
#' @param brlen is the threshold used to distinguish a broad, region-level 
#' events, given in units of fraction of genomic region
regionCall <- function(
	x, 
	brlen = 0.7
){
	x <- lapply(X = x, FUN = "*", attr(x, "pct.overlap"))
	rc_amp <- rc_del <- matrix(
		data = 0L, 
		nrow = nrow(x$pct.amp), 
		ncol = ncol(x$pct.amp), 
		dimnames = dimnames(x$pct.amp)
	)
	rc_amp[is.na(x$pct.amp)] <- NA_integer_
	rc_del[is.na(x$pct.del)] <- NA_integer_
	rc_amp[x$pct.amp > brlen] <-  1L
	rc_amp[x$pct.del > brlen] <- -1L
	rc_del[x$pct.del > brlen] <- -1L
	rc_del[x$pct.amp > brlen] <-  1L
	if(!all(rc_del == rc_amp, na.rm = TRUE)) stop(
		"regionCall: binary status detected for >1 of the regions. ", 
		"Consider increasing 'brlen'"
	)
	return(rc_amp)
}


# Additional support for QDNAseqSignals objects
suppressPackageStartupMessages(require(SummarizedExperiment))
suppressPackageStartupMessages(require(QDNAseq))
#' Get an assay from an eSet object
setMethod(
	f = "assay", 
	signature = c("eSet", "missing"), 
	definition = function(x, i, withDimnames = TRUE, ...){
		assays <- Biobase::assayData(object)
		if(0L == length(assays)) stop(
			"'assay(<", class(x), ">, i=\"missing\", ...) ", 
			"length(assays(<", class(x), ">)) is 0'"
		)
		assays[[1]]
	}
)
setMethod(
	f = "assay", 
	signature = c("eSet", "numeric"), 
	definition = function(x, i, withDimnames = TRUE, ...){
		msg <- paste0(
			"'assay(<", class(x), ">, i=\"numeric\", ...)' ", 
			"invalid subscript 'i'"
		)
		tryCatch(
			expr = {
				Biobase::assayData(object)[[i]]
			}, 
			error = function(err){
				stop(msg, "\n", conditionMessage(err))
			}
		)
	}
)
setMethod(
	f = "assay", 
	signature = c("eSet", "character"), 
	definition = function(x, i, withDimnames = TRUE, ...){
		assays <- Biobase::assayData(object)
		msg <- paste0(
			"'assay(<", class(x), ">, i=\"character\", ...)' ", 
			"invalid subscript 'i'"
		)
		res <- tryCatch(
			expr = {
				Biobase::assayData(object)[[i]]
			}, 
			error = function(err){
				stop(msg, "\n", conditionMessage(err))
			}
		)
		if(is.null(res)) stop(
			msg, "\n'", i, "' not in names(assays(<", class(x), ">))"
		)
		res
	}
)
#' Get the rowRanges from a QDNAseqSignals object
setMethod(
	f = "rowRanges", 
	signature = c("QDNAseqSignals"), 
	definition = function(x, ...){
		as(Biobase::pData(Biobase::featureData(object)), "GRanges")
	}
)
