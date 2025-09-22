# Script for regulatory.MotifAffinityTest
{% include biopipen_dir + "/scripts/regulatory/motifs-common.R" %}

library(BiocParallel)
library(BSgenome)
library(biopipen.utils)

motiffile <- {{in.motiffile | r}}
varfile <- {{in.varfile | r}}
outdir <- {{out.outdir | r}}
ncores <- {{envs.ncores | r}}
tool <- {{envs.tool | r}}
bcftools <- {{envs.bcftools | r}}
genome <- {{envs.genome | r}}
motif_col <- {{envs.motif_col | r}}
regulator_col <- {{envs.regulator_col | r}}
var_col <- {{envs.var_col | r}}
notfound <- {{envs.notfound | r}}
motifdb <- {{envs.motifdb | r}}
regmotifs <- {{envs.regmotifs | r}}
devpars <- {{envs.devpars | r}}
plot_nvars <- {{envs.plot_nvars | r}}
plots <- {{envs.plots | r}}
cutoff <- {{envs.cutoff | r}}
set.seed(8525)

if (is.null(motifdb) || !file.exists(motifdb)) {
    stop("Motif database (envs.motifdb) is required and must exist")
}

if (is.null(genome)) {
    stop("Reference genome (envs.ref) is required and must exist")
}

if (is.null(motiffile) || !file.exists(motiffile)) {
    stop("Motif file (in.motiffile) is required and must exist")
}

if (is.null(varfile) || !file.exists(varfile)) {
    stop("Variant file (in.varfile) is required and must exist")
}

if (is.null(motif_col) && is.null(regulator_col)) {
    stop("Either motif (envs.motif_col) or regulator (envs.regulator_col) column must be provided")
}

log <- get_logger()

log$info("Reading input regulator/motif file ...")
in_motifs <- read.table(motiffile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)


log$info("Ensuring motifs and regulators in the input data ...")
in_motifs <- ensure_regulator_motifs(in_motifs, outdir, motif_col, regulator_col, var_col, regmotifs, notfound = notfound)
genome_pkg <- get_genome_pkg(genome)

motif_var_pairs <- NULL
if (!is.null(var_col)) {
    log$info("Obtaining motif-variant pairs to test ...")
    if (!var_col %in% colnames(in_motifs)) {
        stop("Variant column (envs.var_col) not found in the input motif file")
    }

    motif_var_pairs <- unique(paste0(in_motifs[[motif_col]], " // ", in_motifs[[var_col]]))
}

log$info("Reading variant file ...")
if (grepl("\\.vcf$", varfile) || grepl("\\.vcf\\.gz$", varfile)) {
    log$info("Converting VCF file to BED file ...")
    varfile_bed <- file.path(outdir, gsub("\\.vcf(\\.gz)?$", ".bed", basename(varfile)))
    cmd <- c(
        bcftools, "query",
        "-f", "%CHROM\\t%POS0\\t%END\\t%ID\\t0\\t+\\t%REF\\t%ALT{0}\\n",
        "-i", 'FILTER="PASS" || FILTER="." || FILTER=""',
        "-o", varfile_bed,
        varfile
    )
    run_command(cmd, fg = TRUE)

    varfile <- varfile_bed
}

# `chrom`, `start`, `end`, `name`, `score`, `strand`, `ref`, `alt`.
snpinfo <- read.table(varfile, header=FALSE, stringsAsFactors=FALSE)
colnames(snpinfo) <- c("chrom", "start", "end", "name", "score", "strand", "ref", "alt")

log$info("Reading motif database ...")
mdb <- read_meme_to_motifdb(motifdb, in_motifs, motif_col, regulator_col, notfound, outdir)

tool <- tolower(tool)
tool <- match.arg(tool, c("motifbreakr", "atsnp"))

{% if envs.tool == "motifbreakr" %}
    motifbreakr_args <- {{envs.motifbreakr_args | r}}
    {% include biopipen_dir + "/scripts/regulatory/MotifAffinityTest_MotifBreakR.R" %}
{% else %}
    atsnp_args <- list_update(
        list(padj_cutoff = TRUE, padj = "BH", p = "Pval_diff"),
        {{envs.atsnp_args | r}}
    )
    {% include biopipen_dir + "/scripts/regulatory/MotifAffinityTest_AtSNP.R" %}
{% endif %}
