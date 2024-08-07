# Script for regulatory.MotifAffinityTest
{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "scripts", "regulatory", "motifs-common.R" | source_r }}

library(BiocParallel)
library(BSgenome)

motiffile <- {{in.motiffile | r}}
varfile <- {{in.varfile | r}}
outdir <- {{out.outdir | r}}
ncores <- {{envs.ncores | r}}
tool <- {{envs.tool | r}}
bcftools <- {{envs.bcftools | r}}
genome <- {{envs.genome | r}}
motif_col <- {{envs.motif_col | r}}
regulator_col <- {{envs.regulator_col | r}}
notfound <- {{envs.notfound | r}}
motifdb <- {{envs.motifdb | r}}
regmotifs <- {{envs.regmotifs | r}}
devpars <- {{envs.devpars | r}}
plot_nvars <- {{envs.plot_nvars | r}}
plots <- {{envs.plots | r}}
cutoff <- {{envs.cutoff | r}}

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

log_info("Reading input regulator/motif file ...")
in_motifs <- read.table(motiffile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

log_info("Ensuring motifs and regulators in the input data ...")
in_motifs <- ensure_regulator_motifs(in_motifs, outdir, motif_col, regulator_col, regmotifs, notfound = notfound)
genome_pkg <- get_genome_pkg(genome)

log_info("Reading variant file ...")
if (grepl("\\.vcf$", varfile) || grepl("\\.vcf\\.gz$", varfile)) {
    log_info("Converting VCF file to BED file ...")
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

log_info("Reading motif database ...")
mdb <- read_meme_to_motifdb(motifdb, in_motifs, motif_col, regulator_col, notfound, outdir)

tool <- tolower(tool)
tool <- match.arg(tool, c("motifbreakr", "atsnp"))

if (tool == "motifbreakr") {
    motifbreakr_args <- {{envs.motifbreakr_args | r}}
    {{ biopipen_dir | joinpaths: "scripts", "regulatory", "MotifAffinityTest_MotifBreakR.R" | source_r }}
} else {  # atsnp
    atsnp_args <- {{envs.atsnp_args | r}}
    {{ biopipen_dir | joinpaths: "scripts", "regulatory", "MotifAffinityTest_AtSNP.R" | source_r }}
}
