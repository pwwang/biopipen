library(rlang)
library(universalmotif)
library(MotifDb)
library(biopipen.utils)

#' @title Common functions for regulatory analysis
#' @name regulatory-common
#' @author Panwen Wang

#' Read a regulator-motif mapping file
#'
#' @param rmfile Regulator-motif mapping file
#' @param motif_cols_allowed Allowed motif columns
#' @param reg_cols_allowed Allowed regulator columns
#' @return Data frame with regulators and motifs in the first and second columns, respectively
.read_regmotifs <- function(
    rmfile,
    motif_cols_allowed = c("Motif", "motif", "MOTIF", "Model", "model", "MODEL"),
    reg_cols_allowed = c("Regulator", "regulator", "REGULATOR", "TF", "tf", "TF")
) {
    if (!file.exists(rmfile)) {
        stop("Regulator-motif mapping file does not exist.")
    }
    regmotifs <- read.table(rmfile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
    rm_motif_col <- intersect(motif_cols_allowed, colnames(regmotifs))
    rm_reg_col <- intersect(reg_cols_allowed, colnames(regmotifs))
    if (length(rm_motif_col) == 0) {
        stop(paste0("No motif column found in the regulator-motif mapping file, provide one of: ", paste(motif_cols_allowed, collapse = ", ")))
    }
    if (length(rm_reg_col) == 0) {
        stop(paste0("No regulator column found in the regulator-motif mapping file, provide one of: ", paste(reg_cols_allowed, collapse = ", ")))
    }
    if (length(rm_motif_col) > 1) {
        stop(paste0("Multiple motif columns found (", paste(rm_motif_col, collapse = ", "), ") in the regulator-motif mapping file, provide only one"))
    }
    if (length(rm_reg_col) > 1) {
        stop(paste0("Multiple regulator columns found (", paste(rm_reg_col, collapse = ", "), ") in the regulator-motif mapping file, provide only one"))
    }
    rm_motif_col <- rm_motif_col[1]
    rm_reg_col <- rm_reg_col[1]
    regmotifs <- regmotifs[, c(rm_motif_col, rm_reg_col), drop = FALSE]

    return(regmotifs)
}

#' Handle not found items
#'
#' @param notfound_items Items that were not found
#' @param log_warn Function to log warnings
#' @param msg Message to display
#' @param notfound Action to take if items are not found
#' @param notfound_file File to save the full list of not found items
#' @param log_indent Indentation for log messages
.handle_notfound_items <- function (notfound_items, log_warn, msg, notfound, notfound_file, log_indent = "") {
    if (length(notfound_items) > 0) {
        first_notfound <- head(notfound_items, 3)
        if (length(notfound_items) > 3) {
            first_notfound <- c(first_notfound, "...")
            writeLines(notfound_items, notfound_file)
            msg1 <- paste0(log_indent, msg, ": ", paste(first_notfound, collapse = ", "))
            msg2 <- paste0(log_indent, "Check the full list in ", notfound_file)
            if (notfound == "error") {
                stop(msg1, "\n", msg2)
            } else if (notfound == "ignore") {
                log_warn(msg1)
                log_warn(msg2)
            }
        } else {
            msg <- paste0(log_indent, msg, ": ", paste(first_notfound, collapse = ", "))
            if (notfound == "error") {
                stop(msg)
            } else if (notfound == "ignore") {
                log_warn(msg)
            }
        }
    }
}

#' Read a MEME file to a MotifDb object
#' and filter the motifs based on the input data
#' and return the filtered MotifDb object
#' with metadata
#'
#' @param motifdb MEME file
#' @param indata Input data frame
#' @param motif_col Column name for the motif
#' @param regulator_col Column name for the regulator
#' @param notfound Action to take if motifs are not found
#' @param outdir Output directory, used to save un-matched motifs
#' @return MotifDb object
#' @export
read_meme_to_motifdb <- function(motifdb, indata, motif_col, regulator_col, notfound, outdir) {
    meme <- read_meme(motifdb)
    motifdb_names <- sapply(meme, function(m) m@name)
    motifs <- check_motifs(indata[[motif_col]], motifdb_names, notfound, outdir)
    meme <- filter_motifs(meme, name = motifs)
    # Get the right order of motif names
    motifs <- sapply(meme, function(m) m@name)
    motifdb_matrices <- lapply(meme, function(m) m@motif)
    names(motifdb_matrices) <- motifs
    motifdb_meta <- do.call(rbind, lapply(meme, function(m) {
            ats <- attributes(m)
            ats$dataSource <- strsplit(basename(motifdb), "\\.")[[1]][1]
            ats$class <- NULL
            ats$motif <- NULL
            ats$gapinfo <- NULL
            ats$sequenceCount <- ats$nsites
            ats$providerId <- ats$name
            ats$providerName <- ats$name
            ats$organism <- if (is.null(ats$organism) || length(ats$organism) == 0) "Unknown" else ats$organism
            if (!is.null(regulator_col)) {
                ats$geneSymbol <- indata[
                    indata[[motif_col]] == ats$name,
                    regulator_col,
                    drop = TRUE
                ]
            }
            unlist(ats)
        })
    )
    rownames(motifdb_meta) <- motifs
    MotifDb:::MotifList(motifdb_matrices, tbl.metadata = motifdb_meta)
}

#' Convert a MotifDb object to a motif library
#' with motif names as keys
#' and PWMs as values
#' @param motifdb MotifDb object
#' @return Motif library
#' @export
motifdb_to_motiflib <- function(motifdb) {
    lapply(motifdb, t)
}

#' Make sure the regulators and motifs in the input data from a regulator-motif mappings
#'
#' @param indata Input data frame
#' @param outdir Output directory, used to save un-matched regulators
#' @param motif_col Column name for the motif
#' @param regulator_col Column name for the regulator
#' @param var_col Column name for the variant
#' @param regmotifs Regulator-motif mapping file
#' @param log_indent Indentation for log messages
#' @param notfound Action to take if regulators are not found in the mapping file
#' @return Data frame with regulators and motifs
#' @export
ensure_regulator_motifs <- function (indata, outdir, motif_col, regulator_col, var_col, regmotifs, log_indent = "", notfound = "error", log = NULL) {
    if (is.null(motif_col)) {
        if (is.null(regmotifs)) {
            stop("Regulator-motif mapping file (envs.regmotifs) is required when no motif column (envs.motif_col) is provided")
        }
        log <- log %||% get_logger()
        regmotifs <- .read_regmotifs(regmotifs)
        rm_motif_col <- colnames(regmotifs)[1]
        rm_reg_col <- colnames(regmotifs)[2]
        # check regulators
        rm_regs <- regmotifs[[rm_reg_col]]
        regulators <- indata[[regulator_col]]
        notfound_regs <- setdiff(regulators, rm_regs)
        .handle_notfound_items(
            notfound_regs,
            log$warn,
            "The following regulators were not found in the regulator-motif mapping file",
            notfound,
            file.path(outdir, "notfound_regulators.txt"),
            log_indent
        )
        indata <- indata[indata[[regulator_col]] %in% rm_regs, , drop = FALSE]
        # add motif column
        indata <- merge(indata, regmotifs, by.x = regulator_col, by.y = rm_reg_col, all.x = TRUE, suffixes = c("", "_db"))
        # update motif column
        motif_col <<- rm_motif_col
    } else if (is.null(regulator_col)) {
        if (is.null(regmotifs) || (is.character(regmotifs) && nchar(regmotifs) == 0)) {
            # make motifs unique
            indata <- indata[!duplicated(indata[[motif_col]]), , drop = FALSE]
        } else if (!file.exists(regmotifs)) {
            stop("Regulator-motif mapping file (envs.regmotifs) does not exist.")
        } else {
            # map the regulators
            regmotifs <- .read_regmotifs(regmotifs)
            rm_motif_col <- colnames(regmotifs)[1]
            rm_reg_col <- colnames(regmotifs)[2]
            rm_motifs <- regmotifs[[rm_motif_col]]
            motifs <- indata[[motif_col]]
            notfound_motifs <- setdiff(motifs, rm_motifs)
            .handle_notfound_items(
                notfound_motifs,
                log$warn,
                "The following motifs were not found in the regulator-motif mapping file",
                notfound,
                file.path(outdir, "notfound_motifs.txt"),
                log_indent
            )
            indata <- indata[indata[[motif_col]] %in% rm_motifs, , drop = FALSE]
            # add regulator column
            indata <- merge(indata, regmotifs, by.x = motif_col, by.y = rm_motif_col, all.x = TRUE, suffixes = c("", "_db"))
            # update regulator column
            regulator_col <<- rm_reg_col
        }
    } else {
        indata <- indata[!duplicated(indata[, c(regulator_col, motif_col, var_col), drop = FALSE]), , drop = FALSE]
    }

    return(indata)
}

#' Get the genome package name for a given genome
#'
#' @param genome Genome name
#' @return Genome package name
#' @export
get_genome_pkg <- function(genome) {
    if (!grepl(".", genome, fixed = TRUE)) {
        genome_pkg = sprintf("BSgenome.Hsapiens.UCSC.%s", genome)
    } else {
        genome_pkg = genome
    }
    if (!requireNamespace(genome_pkg, quietly = TRUE)) {
        stop(sprintf("Genome package %s is not installed", genome_pkg))
    }

    library(package = genome_pkg, character.only = TRUE)
    return(genome_pkg)
}

#' Check if motifs are in the motif database
#' and return the motifs that are found
#'
#' @param motifs Motifs to check
#' @param all_motifs All motifs in the motif database
#' @param notfound Action to take if motifs are not found
#' @param outdir Output directory, used to save un-matched motifs
#' @return Motifs that are found
#' @export
check_motifs <- function(motifs, all_motifs, notfound, outdir, log = NULL) {
    log <- log %||% get_logger()
    notfound_motifs <- setdiff(motifs, all_motifs)
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
                log$warn(msg1)
                log$warn(msg2)
            }
        } else {
            msg <- paste0("The following motifs were not found in the motif database: ", paste(first_notfound, collapse = ", "))
            if (notfound == "error") {
                stop(msg)
            } else if (notfound == "ignore") {
                log$warn(msg)
            }
        }

        motifs <- setdiff(motifs, notfound_motifs)
    }
    return(motifs)
}

#' Plot a genomic region surrounding a genomic variant, and
#' potentially disrupted motifs.
#'
#' @param results The motifbreakR results.
#'  A GRanges object with the following columns:
#'  - seqnames: Chromosome
#'  - ranges: Start and end positions
#'  - strand: Strand
#'  -------------------
#'  - SNP_id: Variant ID
#'  - REF: Reference allele
#'  - ALT: Alternative allele
#'  - varType: Variant type. By default, "SNV"
#'  - motifPos: Motif positions
#'  - geneSymbol: Gene symbol, if not provided, try to get from the Regulator column
#'  - dataSource: Motif database source
#'  - providerName: Motif name
#'  - providerId: Motif ID
#'  - effect: Effect of the variant. By default, "strong"
#'  - altPos: Alternative allele position. By default, 1
#'  - alleleDiff: Allele difference, default 0, does not affect the plot for SNVs
#'
#'  Attributes:
#'  - genome.package: Genome package name
#'  - motifs: Motif database, in MotifDb::MotifList format
#' @param variant Variant ID to be plotted
#' @param devpars List of device parameters
#'  - res: Resolution, default 100
#'  - width: Width of the plot, default NULL, calculated based on sequence length
#'  - height: Height of the plot, default NULL, calculated based on the number of motifs
#' @param outdir Output directory. Plots will be saved in the sub-directory "<outdir>/plots/"
#' @export
plot_variant_motifs <- function(results, variant, devpars, outdir) {
    plotdir <- file.path(outdir, "plots")
    dir.create(plotdir, showWarnings = FALSE)

    res <- results[results$SNP_id == variant, , drop = FALSE]
    devpars <- devpars %||% list(res = 100, width = NULL, height = NULL)
    if (length(res) == 0) {
        stop(sprintf("Variant %s not found in results", variant))
    }
    devpars$res <- devpars$res %||% 100
    devpars$height <- devpars$height %||% 2.4 * devpars$res + length(res) * 1.2 * devpars$res
    if (is.null(devpars$width)) {
        left <- min(sapply(res$motifPos, `[`, 1))
        right <- max(sapply(res$motifPos, `[`, 2))
        devpars$width <- 1.5 * devpars$res + (right - left) * 0.3 * devpars$res
        devpars$width <- max(devpars$width, 5 * devpars$res)
    }

    plotfile <- file.path(plotdir, sprintf("%s.png", slugify(variant)))
    # fix motifBreakR 2.12 using names to filter in plotMB
    names(res) <- res$SNP_id
    png(plotfile, width = devpars$width, height = devpars$height, res = devpars$res)
    motifbreakR::plotMB(res, variant)
    dev.off()
}
