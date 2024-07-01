# Script for regulatory.MotifAffinityTest

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

if (is.null(motif_col)) {
    log_info("Inferring motifs from regulators ...")
    if (is.null(regmotifs) || !file.exists(regmotifs)) {
        stop("Regulator motifs (envs.regmotifs) is required and must exist when no motif column (envs.motif_col) is provided")
    }
    regmotifs <- read.table(regmotifs, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    rm_motif_col <- c('Motif', 'motif', 'MOTIF', 'Model', 'model', 'MODEL')
    rm_reg_col <- c('Regulator', 'regulator', 'REGULATOR', 'TF', 'tf', 'TF', 'Transcription factor', 'transcription factor', 'Transcription Factor')
    rm_motif_col <- intersect(rm_motif_col, colnames(regmotifs))
    rm_reg_col <- intersect(rm_reg_col, colnames(regmotifs))
    if (length(rm_motif_col) == 0) {
        stop("No motif column found in envs.regmotifs, provide one of: ", paste(rm_motif_col, collapse = ", "))
    }
    if (length(rm_reg_col) == 0) {
        stop("No regulator column found in envs.regmotifs, provide one of: ", paste(rm_reg_col, collapse = ", "))
    }
    rm_motif_col <- rm_motif_col[1]
    rm_reg_col <- rm_reg_col[1]
    # check regulators
    rm_regs <- regmotifs[, rm_reg_col, drop = TRUE]
    regulators <- in_motifs[, regulator_col, drop = TRUE]
    notfound_regs <- setdiff(regulators, rm_regs)
    if (length(notfound_regs) > 0 && notfound == "error") {
        first_notfound <- head(notfound_regs, 3)
        if (length(notfound_regs) > 3) {
            first_notfound <- c(first_notfound, "...")
            notfound_file <- file.path(outdir, "notfound_regulators.txt")
            writeLines(notfound_regs, notfound_file)
            msg1 <- paste0("The following regulators were not found in the envs.regmotifs file: ", paste(first_notfound, collapse = ", "))
            msg2 <- paste0("Check the full list in ", notfound_file)
            stop(msg1, "\n", msg2)
        } else {
            msg <- paste0("The following regulators were not found in the regmotifs file: ", paste(first_notfound, collapse = ", "))
            stop(msg)
        }
    }
    in_motifs <- in_motifs[in_motifs[, regulator_col] %in% rm_regs, , drop = FALSE]
    # add motif column
    in_motifs <- merge(in_motifs, regmotifs, by.x = regulator_col, by.y = rm_reg_col, all.x = TRUE, suffixes = c("", "_db"))
    motif_col <- rm_motif_col
}
if (is.null(regulator_col)) {
    # make motifs unique
    in_moitfs <- in_motifs[!duplicated(in_motifs[, motif_col]), , drop = FALSE]
} else {
    in_motifs <- in_motifs[!duplicated(in_motifs[, c(regulator_col, motif_col)]), , drop = FALSE]
}


if (!grepl(".", genome, fixed = TRUE)) {
    genome_pkg = sprintf("BSgenome.Hsapiens.UCSC.%s", genome)
} else {
    genome_pkg = genome
}
if (!requireNamespace(genome_pkg, quietly = TRUE)) {
    stop(sprintf("Genome package %s is not installed", genome_pkg))
}

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
        results <- results[1:min(plot_nvars, length(results)), , drop = FALSE]
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
    {% set sourcefile = biopipen_dir | joinpaths: "scripts", "regulatory", "MotifAffinityTest_MotifBreakR.R" %}
    # {{ sourcefile | getmtime }}
    source("{{sourcefile}}")
} else {  # atsnp
    atsnp_args <- {{envs.atsnp_args | r}}
    {% set sourcefile = biopipen_dir | joinpaths: "scripts", "regulatory", "MotifAffinityTest_AtSNP.R" %}
    # {{ sourcefile | getmtime }}
    source("{{sourcefile}}")
}
