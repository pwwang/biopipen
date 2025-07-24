suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
	library(plotthis)
	library(biopipen.utils)
})

indir    <- {{in.indir | r}}
outdir   <- {{out.outdir | r}}
plink    <- {{envs.plink | r}}
indep    <- {{envs.indep | r}}
highld   <- {{envs.highld | r}}
devpars  <- {{envs.devpars | r}}
pihat    <- {{envs.pihat | r}}
samid    <- {{envs.samid | r}}
annofile <- {{envs.anno | r}}
doplot   <- {{envs.plot | r}}
seed     <- {{envs.seed | r}}
ncores   <- {{envs.ncores | r}}

log <- get_logger()

bedfile <- Sys.glob(file.path(indir, '*.bed'))
if (length(bedfile) == 0)
    stop("No bed files found in the input directory.")
if (length(bedfile) > 1) {
    log$warn("Multiple bed files found in the input directory. Using the first one.")
    bedfile <- bedfile[1]
}
input  <- tools::file_path_sans_ext(bedfile)
output <- file.path(outdir, basename(input))

cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--indep-pairwise", indep,
	"--keep-allele-order",
	# One should be mindful of running this with < 50 samples
	# "--bad-ld",
    "--out", output
)
if (!is.null(highld) && !isFALSE(highld)) {
    cmd <- c(cmd, "--range", "--exclude", highld)
}
run_command(cmd, fg = TRUE)

prunein <- paste0(output, '.prune.in')
cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--extract", prunein,
	"--keep-allele-order",
    "--genome",
    "--out", output
)
run_command(cmd, fg = TRUE)

genome <- read.table(
    paste0(output, '.genome'),
    row.names = NULL,
    header = TRUE,
    check.names = FALSE
)
# "unmelt" it
# FID1 IID1 FID2 IID2 RT EZ Z0     Z1     Z2     PI_HAT PHE DST      PPC    RATIO
# s1   s1   s2   s2   UN NA 1.0000 0.0000 0.0000 0.0000 -1  0.866584 0.0000 0.9194
# s1   s1   s2   s2   UN NA 0.4846 0.3724 0.1431 0.3293 -1  0.913945 0.7236 2.0375
# s1   s1   s3   s3   UN NA 1.0000 0.0000 0.0000 0.0000 -1  0.867186 0.0000 1.0791
genome$SAMPLE1 <- paste(genome$FID1, genome$IID1, sep = "\t")
genome$SAMPLE2 <- paste(genome$FID2, genome$IID2, sep = "\t")


# get all samples
samples <- unique(c(genome$SAMPLE1, genome$SAMPLE2))
# make paired into a distance-like matrix
similarity <- genome %>%
    select(SAMPLE1, SAMPLE2, PI_HAT) %>%
    pivot_wider(names_from = SAMPLE2, values_from = PI_HAT, values_fill = NA) %>%
    as.data.frame() %>%
    column_to_rownames("SAMPLE1")
rm(genome)
# get the rownames back
samids <- rownames(similarity)
# get samples that didn't involved
missedrow <- setdiff(samples, rownames(similarity))
missedcol <- setdiff(samples, colnames(similarity))
similarity[missedrow, ] <- NA
similarity[, missedcol] <- NA
# order the matrix
similarity <- similarity[samples, samples, drop = FALSE]
# transpose the matrix to get the symmetric values
sim2 <- t(similarity)
isna <- is.na(similarity)
# fill the na's with their symmetric values
similarity[isna] <- sim2[isna]
rm(sim2)
# still missing: keep them
similarity[is.na(similarity)] <- 0
# get the marks (samples that fail the pihat cutoff)
nsams <- length(samples)
fails <- which(similarity > pihat)
marks <- data.frame(x = (fails - 1)%%nsams + 1, y = ceiling(fails/nsams))
diag(similarity) <- 1

failflags <- rep(F, nrow(marks))
freqs <- as.data.frame(table(factor(as.matrix(marks))))
freqs <- freqs[order(freqs$Freq, decreasing = T), 'Var1', drop = T]
ibd.fail <- c()
while (sum(failflags) < nrow(marks)) {
	samidx <- freqs[1]
	ibd.fail <- c(ibd.fail, samples[samidx])
	freqs <- freqs[-1]
	sapply(1:nrow(marks), function(i) {
		if (samidx %in% marks[i,])
			failflags[i] <<- TRUE
	})
}

ibd_fail_file <- paste0(output, '.ibd.fail')
writeLines(ibd.fail, ibd_fail_file)
cmd <- c(
    plink,
    "--threads", ncores,
    "--bfile", input,
    "--remove", ibd_fail_file,
	"--keep-allele-order",
	"--make-bed",
    "--out", output
)
run_command(cmd, fg = TRUE)

if (doplot) {
	set.seed(seed)
	library(ComplexHeatmap)
	fontsize8 <- gpar(fontsize = 8)
	fontsize9 <- gpar(fontsize = 9)
	ht_opt$heatmap_row_names_gp <- fontsize8
	ht_opt$heatmap_column_names_gp <- fontsize8
	ht_opt$legend_title_gp <- fontsize9
	ht_opt$legend_labels_gp <- fontsize8
	ht_opt$simple_anno_size <- unit(3, "mm")

	samids <- sapply(samples, function(sid) {
		fidiid <- unlist(strsplit(sid, "\t", fixed = TRUE))
		gsub(
            "{fid}",
            fidiid[1],
            gsub("{iid}", fidiid[2], samid, fixed = TRUE),
            fixed = TRUE
        )
	})
	rownames(similarity) <- samids
	colnames(similarity) <- samids

	andata <- NULL
	column_annotation <- NULL
	if (!is.null(annofile) && !isFALSE(annofile)) {
		options(stringsAsFactors = TRUE)
		andata <- read.table(annofile, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
		andata <- andata[samids, , drop = FALSE]
		for (anname in colnames(andata)) {
			column_annotation <- c(column_annotation, anname)
		}
	}

	p <- plotthis::Heatmap(
		similarity,
		name = "PI_HAT",
        in_form = "matrix",
        cell_type = "label",
		rows_data = andata,
		label = function(x) ifelse (x > pihat, '*', NA),
        title = paste0("(*) PI_HAT > ", pihat),
		clustering_distance_rows = function(m) as.dist(1-m),
		clustering_distance_columns = function(m) as.dist(1-m),
        show_row_names = TRUE,
        show_column_names = TRUE,
		column_annotation = column_annotation
	)

	res <- 100
	height <- attr(p, "height") * res
	width <- attr(p, "width") * res
	png(
		filename = paste0(output, '.ibd.png'),
		width = width,
		height = height,
		res = res
	)
	print(p)
	dev.off()

}
