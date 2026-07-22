# CellTypeAnnotation-sctype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R
# Depends on sctype.R being source'd first (for gene_sets_prepare / sctype_score)

annotate_sctype <- function(sobj, ident, tissue, db) {
    library(HGNChelper)

    log <- get_logger()

    if (is.null(db)) { stop("`sctype_db` is not set") }

    # prepare gene sets
    log$info("Preparing gene sets...")
    gs_list <- gene_sets_prepare(db, tissue)

    scRNAseqData <- GetAssayData(sobj, layer = "scale.data")
    idents <- as.character(unique(sobj@meta.data[[ident]]))
    idents <- idents[order(as.numeric(idents))]

    log$info("Working on different levels of cell type labels ...")
    cell_types_list <- list()
    for (i in seq_along(gs_list)) {
        log$info("- Working on level {i} ...")
        if (is.null(gs_list[[i]])) next

        log$info("  Calculating cell-type scores ...")
        es.max <- sctype_score(
            scRNAseqData = scRNAseqData,
            scaled = TRUE,
            gs = gs_list[[i]]$gs_positive,
            gs2 = gs_list[[i]]$gs_negative
        )

        log$info("  Merging cell-type scores by cluster ...")
        cl_resutls <- do_call(
            "rbind",
            lapply(
                idents,
                function(cl) {
                    es.max.cl <- sort(
                        rowSums(es.max[
                            ,
                            rownames(sobj@meta.data[
                                sobj@meta.data[[ident]] == cl,
                            ])
                        ]),
                        decreasing = !0
                    )
                    head(data.frame(
                        cluster = cl,
                        type = names(es.max.cl),
                        scores = es.max.cl,
                        ncells = sum(sobj@meta.data[[ident]] == cl)
                    ), 10)
                }
            )
        )

        sctype_scores <- cl_resutls %>%
            group_by(cluster) %>%
            slice_max(scores, n = 1, with_ties = TRUE)

        if (nrow(sctype_scores) > length(idents)) {
            sctype_scores_count <- sctype_scores %>%
                count(cluster) %>%
                filter(n > 1)
            write("\n########## sctype_scores ###########", stderr())
            write(capture.output(sctype_scores), stderr())
            write("\n####### sctype_scores_count ########", stderr())
            write(capture.output(sctype_scores_count), stderr())
            write("\n####################################", stderr())
            log$info(
                "  Scores tied in the above clusters.",
                immediate. = TRUE
            )
        }

        if (length(gs_list) == 1 || i > 1) {
            # set low-confident (low ScType score) clusters to "unknown"
            log$info("  Setting low-confident clusters to 'Unknown'...")
            sctype_scores$type[
                as.numeric(as.character(sctype_scores$scores)) <
                    sctype_scores$ncells / 4
            ] <- "Unknown"
        }

        celltypes <- sapply(
            idents,
            function(cl) {
                cl_type <- sctype_scores[sctype_scores$cluster == cl, ]
                as.character(cl_type$type[1])
            }
        )
        names(celltypes) <- idents
        cell_types_list[[i]] <- celltypes
    }

    if (length(cell_types_list) == 1) {
        celltypes <- as.list(cell_types_list[[1]])
    } else {
        log$info("Merging cell types at all levels ...")
        celltypes <- list()

        for (i in idents) {
            celltypes[[i]] <- ""
            for (j in seq_along(cell_types_list)) {
                idt <- cell_types_list[[j]][[i]]
                if (idt != "Unknown") {
                    celltypes[[i]] <- paste(celltypes[[i]], idt)
                }
            }
        }
    }

    list(mapping = celltypes)
}
