# Script for regulation.MotifAffinityTest

source("{{biopipen_dir}}/utils/misc.R")
library(BiocParallel)
library(BSgenome)
library(universalmotif)

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

if (!grepl(".", genome, fixed = TRUE)) {
    genome_pkg = sprintf("BSgenome.Hsapiens.UCSC.%s", genome)
} else {
    genome_pkg = genome
}
if (!requireNamespace(genome_pkg, quietly = TRUE)) {
    stop(sprintf("Genome package %s is not installed", genome_pkg))
}

log_info("Reading motif file ...")
in_motifs <- read.table(motiffile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

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
meme <- read_meme(motifdb)

check_motifs <- function(motifdb_names) {
    motifs <- in_motifs[, motif_col, drop = TRUE]
    notfound_motifs <- setdiff(motifs, motifdb_names)
    if (length(notfound_motifs) > 0) {
        first_notfound <- head(notfound_motifs, 3)
        if (length(notfound_motifs) > 3) {
            first_notfound <- c(first_notfound, "...")
            notfound_file <- file.path(outdir, "notfound_motifs.txt")
            writeLines(notfound_motifs, notfound_file)
            msg1 <- paste0("The following motifs were not found in the motif database: ", paste(first_notfound, collapse = ", "))
            msg2 <- paste0("Check the full list in ", notfound_file)

            if (notfound == "error") {
                stop(msg1, "\n", msg2)
            } else if (notfound == "ignore") {
                log_warn(msg1)
                log_warn(msg2)
            }
        } else {
            msg <- paste0("The following motifs were not found in the motif database: ", paste(first_notfound, collapse = ", "))
            if (notfound == "error") {
                stop(msg)
            } else if (notfound == "ignore") {
                log_warn(msg)
            }
        }

        motifs <- setdiff(motifs, notfound_motifs)
    }
    return(motifs)
}

plot_variant <- function(motifbreakr_results) {
    log_info("Plotting variants ...")
    plotdir <- file.path(outdir, "plots")
    dir.create(plotdir, showWarnings = FALSE)
    results <- motifbreakr_results
    if (is.null(plots) || length(plots) == 0) {
        results <- results[order(-abs(results$alleleDiff)), , drop = FALSE]
        results <- results[1:plot_nvars, , drop = FALSE]
        variants <- unique(results$SNP_id)
    } else {
        variants <- names(plots)
    }
    for (variant in variants) {
        log_info("- Variant: {variant}")
        if (is.null(plots[[variant]])) {
            plots[[variant]] <- list(devpars = devpars, which = "TRUE")
        }
        if (is.null(plots[[variant]]$which)) {
            plots[[variant]]$which <- "TRUE"
        }
        if (is.null(plots[[variant]]$devpars)) {
            plots[[variant]]$devpars <- devpars
        }
        if (is.null(plots[[variant]]$devpars$res)) {
            plots[[variant]]$devpars$res <- 100
        }
        res <- results[results$SNP_id == variant, , drop = FALSE]
        if (length(res) == 0) {
            stop(sprintf("Variant %s not found in results", variant))
        }
        res <- subset(res, subset = eval(parse(text = plots[[variant]]$which)))
        if (length(res) == 0) {
            stop(sprintf("No variants to plot for %s", variant))
        }
        plotfile <- file.path(plotdir, sprintf("%s.png", slugify(variant)))
        # fix motifBreakR 2.12 using names to filter in plotMB
        names(res) <- res$SNP_id
        dv <- plots[[variant]]$devpars
        if (is.null(dv$height)) {
            dv$height <- 2.4 * dv$res + length(res) * 1.2 * dv$res
        }
        if (is.null(dv$width)) {
            left <- min(sapply(res$motifPos, `[`, 1))
            right <- max(sapply(res$motifPos, `[`, 2))
            dv$width <- 1.5 * dv$res + (right - left) * 0.3 * dv$res
        }
        png(plotfile, width = dv$width, height = dv$height, res = dv$res)
        motifbreakR::plotMB(res, variant)
        dev.off()
    }
}

tool <- tolower(tool)
tool <- match.arg(tool, c("motifbreakr", "atsnp"))

if (tool == "motifbreakr") {
    motifbreakr_args <- {{envs.motifbreakr_args | r}}
    {% set sourcefile = biopipen_dir | joinpaths: "scripts", "regulation", "MotifAffinityTest_MotifBreakR.R" %}
    # {{ sourcefile | getmtime }}
    source("{{sourcefile}}")
} else {  # atsnp
    atsnp_args <- {{envs.atsnp_args | r}}
    {% set sourcefile = biopipen_dir | joinpaths: "scripts", "regulation", "MotifAffinityTest_AtSNP.R" %}
    # {{ sourcefile | getmtime }}
    source("{{sourcefile}}")
}
