# CellTypeAnnotation-garnett.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R and GarnettClassifierTrainer.R

# patch_garnett_make_predictions — replacement for the wrapper function of the
# same name in biopipen/scripts/scrna/CellTypeAnnotation-garnett.R
#
# Root cause being fixed: with glmnet >= 4.0, predict(cv.glmnet, type="response")
# on a multinomial fit returns a rank-3 array [cells, classes, s]. garnett's
# `as.matrix(as.data.frame(temp))` then pastes dim2 and dim3 names together, so
# the class columns come out as "B cell.lambda.min" instead of "B cell", and
# every class name is mangled before the which.max()/gate logic runs. Fix: cut
# the singleton `s` dimension down to a plain cells x classes matrix before any
# data.frame conversion, strip any residual ".lambda.<s>" / ".s=<x>" / ".1"
# suffix, and refuse loudly if the class dimension cannot be established
# (all-Unknown is a legitimate result; silently collapsing the classes is not).
#
# Usage: replace the existing patch_garnett_make_predictions() in
# CellTypeAnnotation-garnett.R with this function. The rest of the wrapper is
# unchanged (it is still called once at the top of annotate_garnett()).
patch_garnett_make_predictions <- function(log) {
    make_predictions_fixed <- function(cds, classifier, curr_node,
                                       rank_prob_ratio, cores = 1, s) {
        cvfit <- igraph::V(classifier@classification_tree)[curr_node]$model[[1]]
        predictions <- tryCatch({
            if (is.null(cvfit)) {
                child_cell_types <- igraph::V(
                    classifier@classification_tree
                )[suppressWarnings(.outnei(curr_node))]$name
                predictions <- matrix(
                    FALSE,
                    nrow = nrow(colData(cds)),
                    ncol = length(child_cell_types),
                    dimnames = list(row.names(colData(cds)), child_cell_types)
                )
                predictions <- split(
                    predictions,
                    rep(1:ncol(predictions), each = nrow(predictions))
                )
                names(predictions) <- child_cell_types
                predictions
            } else {
                candidate_model_genes <-
                    cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
                good_genes <- intersect(
                    row.names(counts(cds)), candidate_model_genes
                )
                if (length(good_genes) == 0) {
                    stop(paste(
                        "None of the model genes are in your CDS object.",
                        "Did you specify the correct cds_gene_id_type and",
                        "the correct db?"
                    ))
                }
                x <- Matrix::t(counts(cds[
                    intersect(
                        row.names(counts(cds)), candidate_model_genes
                    ),
                ]))
                extra <- as(matrix(
                    0,
                    nrow = nrow(x),
                    ncol = length(setdiff(candidate_model_genes, colnames(x)))
                ), "sparseMatrix")
                row.names(extra) <- row.names(x)
                colnames(extra) <- setdiff(candidate_model_genes, colnames(x))
                x <- cbind(x, extra)
                x <- x[, candidate_model_genes]
                nonz <- Matrix::rowSums(do.call(
                    cbind, glmnet::coef.glmnet(cvfit, s = "lambda.min")
                ))
                nonz <- nonz[2:length(nonz)]
                nonz <- names(nonz[nonz != 0])
                if (sum(!nonz %in% row.names(counts(cds))) > 0) {
                    warning(paste(
                        "The following genes used in the classifier are not",
                        "present in the input CDS. Interpret with caution.",
                        nonz[!nonz %in% row.names(counts(cds))]
                    ))
                }
                temp <- stats::predict(cvfit, newx = x, s = s, type = "response")
                temp[is.nan(temp)] <- 0
                # --- FIX (a): with glmnet >= 4 the multinomial response is a
                # [cells, classes, s] array. Rebuild it as a plain cells x
                # classes matrix before any data.frame conversion. Two traps:
                # `temp[, , 1, drop = FALSE]` keeps rank 3, so as.data.frame()
                # below pastes dim2+dim3 names together and every class comes
                # out as "<class>.<s-label>"; plain `temp[, , 1]` instead drops
                # any length-1 dim, so a single-class model collapses to a bare
                # vector and the class name is lost. Build the matrix by hand
                # so exactly the class labels survive either way.
                if (is.array(temp) && length(dim(temp)) > 2) {
                    if (length(dim(temp)) != 3) {
                        stop(paste(
                            "Unexpected predict() rank",
                            length(dim(temp)), "(dim =",
                            paste(dim(temp), collapse = " x "),
                            "); expected [cells, classes, s]"
                        ))
                    }
                    if (dim(temp)[3] != 1) {
                        stop(paste(
                            "Expected a single `s` slot in the predict() output,",
                            "got", dim(temp)[3], "(dim =",
                            paste(dim(temp), collapse = " x "),
                            "); refusing to guess the class dimension"
                        ))
                    }
                    temp <- matrix(
                        temp[, , 1],
                        nrow = dim(temp)[1],
                        ncol = dim(temp)[2],
                        dimnames = list(dimnames(temp)[[1]], dimnames(temp)[[2]])
                    )
                }
                if (is.null(temp) || is.null(dim(temp)) || ncol(temp) < 2) {
                    stop(paste(
                        "The classifier returned no usable class dimension",
                        "(dim =", paste(dim(temp), collapse = " x "),
                        "); cannot classify"
                    ))
                }
                colnames(temp) <- sub("\\.lambda\\..*$", "", colnames(temp))
                colnames(temp) <- sub("\\.s=.*$", "", colnames(temp))
                colnames(temp) <- sub("\\.1$", "", colnames(temp))
                prediction_probs <- as.matrix(as.data.frame(temp))
                prediction_probs <-
                    prediction_probs / Biobase::rowMax(prediction_probs)
                prediction_probs[is.nan(prediction_probs)] <- 0
                prediction_probs <- apply(prediction_probs, 1, function(x) {
                    m <- names(which.max(x))
                    s <- sort(x, decreasing = T)
                    c(cell_type = m, odds_ratio = s[1] / s[2])
                })
                prediction_probs <- as.data.frame(t(prediction_probs))
                prediction_probs$cell_name <- row.names(prediction_probs)
                names(prediction_probs) <- c(
                    "cell_type", "odds_ratio", "cell_name"
                )
                prediction_probs$odds_ratio <- as.numeric(
                    as.character(prediction_probs$odds_ratio)
                )
                assignments <- prediction_probs[
                    prediction_probs$odds_ratio > rank_prob_ratio,
                ]
                random_guess_thresh <- 1 / length(cvfit$glmnet.fit$beta)
                assignments <- assignments[
                    assignments$odds_ratio > random_guess_thresh,
                ]
                not_assigned <- row.names(colData(cds))[
                    !row.names(colData(cds)) %in% assignments$cell_name
                ]
                if (length(not_assigned) > 0) {
                    assignments <- rbind(
                        assignments,
                        data.frame(
                            cell_name = not_assigned,
                            cell_type = NA, odds_ratio = NA
                        )
                    )
                }
                assignments$cell_type <- stringr::str_replace_all(
                    assignments$cell_type, "\\.1", ""
                )
                predictions <- reshape2::dcast(
                    assignments,
                    cell_name ~ cell_type,
                    value.var = "odds_ratio"
                )
                predictions <- predictions[!is.na(predictions$cell_name), ]
                row.names(predictions) <- predictions$cell_name
                if (ncol(predictions) > 2) {
                    predictions <- predictions[
                        ,
                        setdiff(colnames(predictions), "NA")
                    ]
                    predictions <- predictions[, -1, drop = FALSE]
                    predictions <- predictions[
                        rownames(colData(cds)), , drop = FALSE
                    ]
                    predictions <- as.matrix(predictions)
                    predictions[is.na(predictions)] <- FALSE
                    predictions[predictions != 0] <- TRUE
                    cell_type_names <- colnames(predictions)
                    predictions <- split(
                        predictions,
                        rep(1:ncol(predictions), each = nrow(predictions))
                    )
                    names(predictions) <- cell_type_names
                } else {
                    cell_type_names <- names(cvfit$glmnet.fit$beta)
                    one_type <- names(predictions)[2]
                    if (one_type == "NA") {
                        names(predictions)[2] <- "Unknown"
                        one_type <- "Unknown"
                    }
                    predictions <- matrix(
                        FALSE,
                        nrow = nrow(colData(cds)),
                        ncol = length(cell_type_names),
                        dimnames = list(
                            row.names(colData(cds)), cell_type_names
                        )
                    )
                    # --- FIX (b): "Unknown" (and any name outside the model's
                    # classes) must not index the matrix; leave the mask
                    # all-FALSE instead of crashing.
                    if (one_type %in% cell_type_names) {
                        predictions[, one_type] <- TRUE
                    }
                    predictions <- split(
                        predictions,
                        rep(1:ncol(predictions), each = nrow(predictions))
                    )
                    names(predictions) <- cell_type_names
                }
                predictions
            }
        }, error = function(e) {
            if (e$message == paste(
                "None of the model genes are in your CDS object.",
                "Did you specify the correct cds_gene_id_type and",
                "the correct db?"
            )) {
                stop(e)
            }
            print(e)
            cell_type_names <- names(cvfit$glmnet.fit$beta)
            predictions <- matrix(
                FALSE,
                nrow = nrow(colData(cds)),
                ncol = length(cell_type_names),
                dimnames = list(row.names(colData(cds)), cell_type_names)
            )
            predictions <- split(
                predictions,
                rep(1:ncol(predictions), each = nrow(predictions))
            )
            names(predictions) <- cell_type_names
            predictions
        })
    }
    tryCatch({
        monkey_patch("garnett", "make_predictions", make_predictions_fixed)
    }, error = function(e) {
        log$warn(paste(
            "Failed to patch garnett::make_predictions (glmnet >= 4.0 fix):",
            conditionMessage(e)
        ))
    })
}

# Workaround for garnett 0.2.22: the marker-file lexer (`t_NAME` in the rly
# `Lexer` in garnett's namespace) only accepts ASCII letters, so a cell type
# like "γδ-T cells" from ScTypeDB fails with `Marker file error. Syntax error
# 'γ' ...` when train_cell_classifier() parses the file. Widen the token
# regex to Unicode letters via (*UCP) (ASCII-only behavior is unchanged, as
# [:alpha:]/[:alnum:] still cover the ASCII range) so such names pass through
# verbatim. Idempotent — safe to call before train_cell_classifier().
patch_garnett_marker_lexer <- function(log) {
    tryCatch({
        # parse_input() rebuilds the lexer from the generator at every call
        # (rly::lex(Lexer)), so patching the generator in garnett's namespace
        # is enough; its environment is not locked (R6 lock_class = FALSE)
        ns <- asNamespace("garnett")
        Lexer <- get("Lexer", envir = ns)
        fields <- Lexer$public_fields
        fields[["t_NAME"]] <- "(*UCP)[[:alnum:]_+/\\-\\.|=`~\\*&<^%?@!$();:]*[[:alpha:]][[:alnum:]_+/\\-\\.|=`~\\*&<^%?@!$();]*"
        Lexer$public_fields <- fields
    }, error = function(e) {
        log$warn(paste(
            "Failed to patch garnett's marker-file lexer for Unicode cell",
            "type names:",
            conditionMessage(e)
        ))
    })
}

# patch_garnett_run_classifier — fix the upstream crash in garnett 0.2.22's
# run_classifier() when NO cell passes any gate.
#
# Upstream bug: run_classifier() builds its level table with
#   level_table <- data.frame(cell = row.names(imputed_gate_res[[1]]), ...)
# but make_predictions() returns split()-ed plain vectors, whose row.names() is
# NULL, and data.frame(cell = NULL, level1 = "Unknown") errors in R >= 4.x with
# "arguments imply differing number of rows: 0, 1". The line is unconditional,
# so this guard is what makes a degenerate classifier fail loudly as
# all-Unknown instead of killing the job.
# Fix: take cell names from the cds instead of from the prediction vectors, and
# warn when no cell was classified (degenerate classifier / ratio too strict).
patch_garnett_run_classifier <- function(log) {
    fix_run_classifier <- function() {
        src <- deparse(garnett:::run_classifier)
        n_cells <- "nrow(SummarizedExperiment::colData(cds))"
        cell_names <- "row.names(SummarizedExperiment::colData(cds))"
        src <- gsub(
            "length(imputed_gate_res[[1]])", n_cells, src, fixed = TRUE
        )
        src <- gsub(
            "row.names(imputed_gate_res[[1]])", cell_names, src, fixed = TRUE
        )
        warn_at <- grep("tree_levels <- igraph::distances", src, fixed = TRUE)
        warn_line <- paste0(
            "    if (length(imputed_gate_res) == 0 || ",
            "!any(vapply(imputed_gate_res, length, integer(1)) > 0)) ",
            "warning(\"garnett run_classifier: no cell passed any gate; all ",
            "cells marked Unknown (degenerate classifier or rank_prob_ratio ",
            "too strict)\")"
        )
        if (length(warn_at) > 0) {
            src <- append(src, warn_line, after = warn_at[1] - 1)
        }
        # A silent no-op would be worse than a loud failure: if garnett's
        # internals change, the substitutions above stop matching.
        if (!any(grepl(n_cells, src, fixed = TRUE)) ||
                any(grepl("imputed_gate_res[[1]]", src, fixed = TRUE))) {
            stop(paste(
                "run_classifier substitution did not match the installed",
                "garnett::run_classifier source; the empty-gate crash fix is",
                "NOT applied"
            ))
        }
        eval(
            parse(text = paste(src, collapse = "\n")),
            envir = asNamespace("garnett")
        )
    }
    tryCatch(
        {
            monkey_patch("garnett", "run_classifier", fix_run_classifier())
            log$info("Patched garnett::run_classifier (empty-gate crash fix)")
        },
        error = function(e) {
            log$warn(paste(
                "Failed to patch garnett::run_classifier",
                "(empty-gate crash fix):", conditionMessage(e)
            ))
        }
    )
}

annotate_garnett <- function(sobj, ident, garnett_args) {
    library(monocle3)
    library(garnett)
    library(SeuratWrappers)

    log <- get_logger()
    patch_garnett_make_predictions(log)
    patch_garnett_run_classifier(log)
    patch_garnett_marker_lexer(log)

    classifier_path <- garnett_args$classifier
    if (is.null(classifier_path)) { stop("`envs.garnett.classifier` is not set") }
    if (!file.exists(classifier_path)) {
        stop(paste0("Garnett classifier file does not exist: ", classifier_path))
    }

    log$info("Loading Garnett classifier ...")
    classifier <- read_obj(classifier_path)
    if (!inherits(classifier, "garnett_classifier")) {
        stop(paste0(
            "The file in `envs.garnett.classifier` does not contain a ",
            "Garnett classifier (garnett_classifier) object: ", classifier_path
        ))
    }

    # Gene ID conversion database: "none" or an AnnotationDb/OrgDb package name
    db <- garnett_args$db %||% "none"
    if (!identical(db, "none")) {
        if (!requireNamespace(db, quietly = TRUE)) {
            stop(paste0(
                "The gene ID database package `", db,
                "` from `envs.garnett.db` is not installed."
            ))
        }
        library(db, character.only = TRUE)
        db <- getExportedValue(db, db)
    }

    assay <- garnett_args$assay
    # Keys consumed here must not be forwarded to classify_cells()
    garnett_args$classifier <- NULL
    garnett_args$db <- NULL
    garnett_args$assay <- NULL

    unknown_args <- setdiff(names(garnett_args), formalArgs(classify_cells))
    if (length(unknown_args) > 0) {
        stop(paste0(
            "Unknown arguments in `envs.garnett`: ",
            paste(unknown_args, collapse = ", "),
            " (not arguments of `garnett::classify_cells()`)"
        ))
    }

    log$info("Converting Seurat object to a monocle3 cell_data_set ...")
    # as.cell_data_set is the monocle3 generic; SeuratWrappers only registers
    # the method for Seurat objects (loaded above), so call it unqualified
    cds <- if (is.null(assay)) {
        as.cell_data_set(sobj)
    } else {
        as.cell_data_set(sobj, assay = assay)
    }

    # classify_cells() hard-asserts a Size_Factor column but overwrites the
    # value internally with colSums / (cell_totals * median(num_genes_expressed))
    # -- the monocle3 convention (colSums / G) -- so the labels do not depend on
    # what we pass here; the returned cds does keep it though. Install monocle3
    # (median-ratio) size factors unconditionally so the annotation cds stays on
    # the same convention the trainer uses, instead of the
    # SeuratWrappers::as.cell_data_set() colSums convention, and verify the
    # conversion actually happened (a no-op monocle3 would leave the colSums
    # values in place, which is the silent all-"Unknown" trap at training time).
    cds <- monocle3::estimate_size_factors(cds)
    .sf <- as.numeric(colData(cds)$Size_Factor)
    .ct <- as.numeric(Matrix::colSums(counts(cds)))
    .g <- exp(mean(log(.ct[.ct > 0])))
    if (anyNA(.sf) || !all(is.finite(.sf))) {
        stop("Could not install finite size factors on the cell_data_set.")
    }
    if (!isTRUE(all.equal(.g, 1)) &&
        isTRUE(all.equal(.sf, .ct, check.attributes = FALSE))) {
        stop(paste(
            "Garnett expects monocle3 (median-ratio) size factors, but",
            "colData(cds)$Size_Factor still equals colSums(counts) (the",
            "SeuratWrappers::as.cell_data_set() convention). At training time",
            "this makes the model fit and the prediction matrix differ by a",
            "factor of G = exp(mean(log(colSums))), which turns every cell into",
            "'Unknown' with no error. Run monocle3::estimate_size_factors(cds)",
            "on the converted cell_data_set."
        ))
    }

    # classify_cells() internally normalizes counts(cds); an assay with only a
    # data layer (converted counts are log values) would silently give garbage
    if (is(counts(cds), "dgCMatrix") && any(counts(cds)@x %% 1 != 0)) {
        log$warn(paste(
            "The counts of the converted cell_data_set are not integers;",
            "garnett::classify_cells() expects raw counts. Results may be",
            "unreliable. Set `envs.assay` to an assay with a counts layer."
        ))
    }

    log$info("Classifying cells with Garnett ...")
    classify_args <- garnett_args
    classify_args$cds <- cds
    classify_args$classifier <- classifier
    classify_args$db <- db
    result_cds <- do_call(classify_cells, classify_args)

    labels <- colData(result_cds)$cell_type
    if (is.null(labels)) {
        stop("classify_cells() did not return a `cell_type` column")
    }
    # Cells with zero counts are excluded by classify_cells and get NA
    labels[is.na(labels)] <- "Unknown"

    n_unknown <- sum(labels == "Unknown")
    if (n_unknown == length(labels)) {
        log$warn(paste(
            "All cells are classified as 'Unknown' by Garnett; check",
            "`envs.garnett.db`/`envs.garnett.cds_gene_id_type` (classifier is",
            "trained on", classifier@gene_id_type, "genes)"
        ))
    } else {
        log$info(
            "Garnett classified {length(labels) - n_unknown}/{length(labels)}",
            " cells (Unknown: {n_unknown})"
        )
    }

    result <- data.frame(
        garnett_celltype = unname(labels),
        row.names = colnames(result_cds)
    )

    if (is.null(ident)) {
        list(cell_annotations = result, annotation_col = "garnett_celltype")
    } else {
        log$info("Aggregating Garnett results by cluster ...")
        mapping <- majority_vote(
            labels, as.character(sobj@meta.data[[ident]]),
            unknown = "Unknown"  # garnett's sentinel is capital-U
        )
        list(
            mapping = mapping,
            cell_annotations = result,
            annotation_col = "garnett_celltype"
        )
    }
}
