library(rlang)
library(rtracklayer)
library(MatrixEQTL)
library(biopipen.utils)

snpfile = {{in.geno | r}}
expfile = {{in.expr | r}}
covfile = {{in.cov | r}}
joboutdir = {{job.outdir | r}}
alleqtl = {{out.alleqtls | r}}
outfile = {{out.cisqtls | r}}

model = {{envs.model | r}}
pval = {{envs.pval | r}}
match_samples = {{envs.match_samples | r}}
transp = {{envs.transp | r}}
fdr = {{envs.fdr | r}}
snppos = {{envs.snppos | r}}
genepos = {{envs.genepos | r}}
dist = {{envs.dist | r}}

transpose_geno = {{envs.transpose_geno | r}}
transpose_expr = {{envs.transpose_expr | r}}
transpose_cov = {{envs.transpose_cov | r}}

log <- get_logger()

arg_match(model, c("modelANOVA", "modelLINEAR", "linear", "anova"))
if (model == "linear") model = "modelLINEAR"
if (model == "anova") model = "modelANOVA"
model = get(model)

trans_enabled = !is.null(transp)
cis_enabled = !is.null(snppos) && !is.null(genepos) && dist > 0

# if trans is disabled, all files needed for cis should be provided
if (!trans_enabled && !cis_enabled) {
    log$warn("Using `envs.transp = 1e-5` since cis-eQTL is disabled.")
    trans_enabled <- TRUE
    transp <- 1e-5
}

transpose_file <- function(file, what) {
    if (is.null(file)) return(NULL)
    log$info("Transposing {what} file ...")
    out <- file.path(joboutdir, paste0(
        tools::file_path_sans_ext(basename(file)),
        ".transposed.",
        tools::file_ext(file))
    )
    data <- read.table(file, header=TRUE, stringsAsFactors=FALSE, row.names=1, sep="\t", quote="", check.names=FALSE)
    write.table(t(data), file=out, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
    out
}

if (transpose_geno) snpfile = transpose_file(snpfile, "geno")
if (transpose_expr) expfile = transpose_file(expfile, "expr")
if (transpose_cov) covfile = transpose_file(covfile, "cov")

log$info("Loading SNP data ...")
snps = SlicedData$new();
snps$fileDelimiter = "\t";       # the TAB character
snps$fileOmitCharacters = "NA";  # denote missing values;
snps$fileSkipRows = 1;           # one row of column labels
snps$fileSkipColumns = 1;        # one column of row labels
snps$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps$LoadFile( snpfile );

log$info("Loading gene expression data ...")
gene = SlicedData$new();
gene$fileDelimiter = "\t";       # the TAB character
gene$fileOmitCharacters = "NA";  # denote missing values;
gene$fileSkipRows = 1;           # one row of column labels
gene$fileSkipColumns = 1;        # one column of row labels
gene$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
gene$LoadFile( expfile );

cvrt = SlicedData$new();
if (!is.null(covfile) && file.exists(covfile)) {
    log$info("Loading covariate data ...")
    covmatrix = read.table(covfile, header=TRUE, stringsAsFactors=FALSE, row.names=1, sep="\t", quote="", check.names=FALSE)
    cvrt$CreateFromMatrix( as.matrix(covmatrix) )
}

log$info("Matching samples ...")
if (match_samples) {
    # let matrixEQTL raise an error if samples do not match
} else {
    n_sample_snps = snps$nCols()
    n_sample_gene = gene$nCols()
    common_samples = intersect(snps$columnNames, gene$columnNames)
    if (!is.null(covfile)) {
        common_samples = intersect(common_samples, cvrt$columnNames)
        n_sample_cov = cvrt$nCols()
        cvrt = cvrt$ColumnSubsample(match(common_samples, cvrt$columnNames))
    }
    snps = snps$ColumnSubsample(match(common_samples, snps$columnNames))
    gene = gene$ColumnSubsample(match(common_samples, gene$columnNames))
    log$info("- Samples used in SNP data: {n_sample_snps} -> {snps$nCols()}")
    log$info("- Samples used in gene expression data: {n_sample_gene} -> {gene$nCols()}")
    if (!is.null(covfile)) {
        log$info("- Samples used in covariate data: {n_sample_cov} -> {cvrt$nCols()}")
    }
}

log$info("Composing engine parameters ...")
engine_params = list()
engine_params$snps = snps
engine_params$gene = gene
engine_params$cvrt = cvrt
engine_params$output_file_name = if(trans_enabled) alleqtl else NULL
engine_params$pvOutputThreshold = if(trans_enabled) min(transp, 1) else 0
engine_params$useModel = model
engine_params$errorCovariance = numeric()
engine_params$verbose = TRUE
engine_params$noFDRsaveMemory = !fdr

noq = function(s) {
    gsub('^\"|\"$', "", s)
}

if (cis_enabled) {
    log$info("Loading SNP positions ...")
    if (endsWith(snppos, ".bed")) {
        snppos_data = read.table(snppos, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        snppos_data = data.frame(
            snp = snppos_data$V4,
            chr = snppos_data$V1,
            pos = snppos_data$V3
        )
    } else if (endsWith(snppos, ".gff") || endsWith(snppos, ".gtf")) {
        snppos_data = import(snppos)
        elem_meta = elementMetadata(snppos_data)
        snppos_data = data.frame(
            snp = elem_meta$snp_id %||% elem_meta$rs_id %||% elem_meta$rs,
            chr = as.character(seqnames(snppos_data)),
            pos = start(snppos_data)
        )
    } else if (endsWith(snppos, ".vcf") || endsWith(snppos, ".vcf.gz")) {
        snppos_data = read.table(
            snppos,
            header=FALSE,
            row.names=NULL,
            stringsAsFactors=FALSE,
            check.names=FALSE
        )
        snppos_data = snppos_data[, c(3, 1, 2)]
        colnames(snppos_data) = c("snp", "chr", "pos")
    } else {
        # snp	chr	pos
        # Snp_01	chr1	721289
        # Snp_02	chr1	752565
        # check if 3rd column of the first line is numeric.
        # if it is, there is no header; otherwise, it is a header.
        header <- is.na(suppressWarnings(as.numeric(strsplit(readLines(snppos, n = 1), "\t")[[1]][3])))

        snppos_data = read.table(
            snppos,
            sep = "\t",
            header = header,
            row.names = NULL,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        colnames(snppos_data) = c("snp", "chr", "pos")
    }

    log$info("Loading gene positions ...")
    if (endsWith(genepos, ".bed")) {
        genepos_data = read.table(genepos, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        genepos_data = data.frame(
            geneid = genepos_data$V4,
            chr = genepos_data$V1,
            s1 = genepos_data$V2,
            s2 = genepos_data$V3
        )
    } else if (endsWith(genepos, ".gff") || endsWith(genepos, ".gtf")) {
        genepos_data = import(genepos)
        elem_meta = elementMetadata(genepos_data)
        genepos_data = data.frame(
            geneid = elem_meta$gene_id %||% elem_meta$gene_name,
            chr = as.character(seqnames(genepos_data)),
            s1 = start(genepos_data),
            s2 = end(genepos_data)
        )
    } else {
        parts <- strsplit(readLines(genepos, n = 1), "\t")[[1]]
        header <- is.na(suppressWarnings(as.numeric(parts[3]))) || is.na(suppressWarnings(as.numeric(parts[4])))
        genepos_data = read.table(
            genepos,
            sep = "\t",
            header = header,
            row.names = NULL,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        colnames(genepos_data) = c("geneid", "chr", "s1", "s2")
    }

    log$info("Running MatrixEQTL with cis-eQTLs enabled ...")
    engine_params$output_file_name.cis = outfile
    engine_params$pvOutputThreshold.cis = min(pval, 1)
    engine_params$cisDist = dist
    engine_params$snpspos = snppos_data
    engine_params$genepos = genepos_data
    do_call(Matrix_eQTL_main, engine_params)
    if (!file.exists(alleqtl)) file.create(alleqtl)
} else {
    log$info("Running MatrixEQTL without cis-eQTLs ...")
    do_call(Matrix_eQTL_engine, engine_params)
    if (!file.exists(outfile)) file.create(outfile)
}

if (pval == 0) {
    if (!file.exists(outfile)) file.create(outfile)
    if (!file.exists(alleqtl)) file.create(alleqtl)
}
