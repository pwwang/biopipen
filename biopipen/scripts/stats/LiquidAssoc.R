library(rlang)
library(dplyr)
library(tidyr)
library(fastLiquidAssociation)
library(biopipen.utils)

infile <- {{in.infile | r}}
covfile <- {{in.covfile | r: quote_none=False | r}}
groupfile <- {{in.groupfile | r}}
fmlfile <- {{in.fmlfile | r}}
outfile <- {{out.outfile | r}}
x <- {{envs.x | r}}
nvec <- {{envs.nvec | r}}
topn <- {{envs.topn | r}}
rvalue <- {{envs.rvalue | r}}
cut <- {{envs.cut | r}}
ncores <- {{envs.ncores | r}}
padj <- {{envs.padj | r}}
transpose_input <- {{envs.transpose_input | r}}
transpose_group <- {{envs.transpose_group | r}}
transpose_cov <- {{envs.transpose_cov | r}}
xyz_names <- {{envs.xyz_names | r}}
if (!is.null(xyz_names) && length(xyz_names) == 1) {
	xyz_names <- trimws(strsplit(xyz_names, ",")[[1]])
}

if (is.null(groupfile) && is.null(nvec)) {
	stop("Must provide either in.groupfile or envs.nvec")
}
if (!is.null(groupfile) && !is.null(nvec)) {
	stop("Must provide either in.groupfile or envs.nvec, not both")
}

log$info("Reading and preparing data ...")
indata <- read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
if (transpose_input) {
	indata <- t(indata)
}
if (!is.null(covfile)) {
	covdata <- read.table(covfile, header = TRUE, sep = "\t", row.names = 1)
	if (transpose_cov) {
		covdata <- t(covdata)
	}
	if (!isTRUE(all.equal(rownames(indata), rownames(covdata)))) {
		stop("Row names of indata and covdata must be identical")
	}
	indata <- indata %>% mutate(across(everything(), function(xx) {
		lm(xx ~ as.matrix(covdata))$residuals
	}))
}

expand_range <- function(range) {
	items <- trimws(strsplit(range, ",|-")[[1]])
	num_items <- as.numeric(items)
	if (anyNA(num_items)) {
		# it's sample names
		return(match(items, colnames(indata)))
	}
	return(num_items)
}

cut <- cut %||% max(ceiling(nrow(indata)/22), 4)
if (!is.null(x)) { x <- expand_range(x) }
if (!is.null(groupfile)) {
	groupdata <- read.table(groupfile, header = TRUE, sep = "\t", row.names = 1)
	if (transpose_group) {
		groupdata <- t(groupdata)
	}
	if (!isTRUE(all.equal(rownames(indata), rownames(groupdata)))) {
		stop("Row names of indata and groupdata must be identical")
	}
	nvec <- (ncol(indata) + 1) : (ncol(indata) + ncol(groupdata))
	indata <- cbind(indata, groupdata)
} else {
	nvec <- expand_range(nvec)
}

log$info("Running fastLiquidAssociation ...")
indata <- as.matrix(indata)
mla <- fastMLA(
	data = indata,
	topn = topn,
	rvalue = rvalue,
	cut = cut,
	threads = ncores,
	nvec = nvec
)

if (nrow(mla) == 0) {
	log$warn("No significant associations found")
	out <- data.frame(
		X12 = character(),
		X21 = character(),
		X3 = character(),
		rhodiff = numeric(),
		`MLA.value` = numeric(),
		estimates = numeric(),
		`san.se` = numeric(),
		wald = numeric(),
		Pval = numeric(),
		model = character()
	)
} else {
	cnm <- mass.CNM(data = indata, GLA.mat = mla, nback = topn)
	out <- cnm$`top p-values` %>%
		dplyr::select(X12 = "X1 or X2", X21 = "X2 or X1", everything(), Pval = "p value")
}

if (!is.null(fmlfile)) {
	fmldata <- read.table(fmlfile, header = FALSE, sep = "\t", row.names = NULL)
	colnames(fmldata) <- c("Z", "X", "Y")
	all_combns <- fmldata %>% unite("XYZ", X, Y, Z, sep = " // ") %>% pull(XYZ)
	out <- out %>%
		unite("XYZ", X12, X21, X3, sep = " // ", remove = FALSE) %>%
		dplyr::filter(XYZ %in% all_combns) %>%
		dplyr::select(-XYZ)
}

if (!is.null(xyz_names)) {
	out <- out %>%
		dplyr::select(
			!!sym(xyz_names[1]) := "X12",
			!!sym(xyz_names[2]) := "X21",
			!!sym(xyz_names[3]) := "X3",
			everything()
		)
}

if (padj != "none") {
	log$info("Calculating adjusted p-values ...")
	out$Padj <- p.adjust(out$Pval, method = padj)
}

log$info("Writing output ...")
write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
