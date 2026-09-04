# CellTypeAnnotation-llmcelltype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_llmcelltype <- function(sobj, ident, llmcelltype_args) {
    library(LLMCellType)

    log <- get_logger()

    llmcelltype_argnames <- formalArgs(LLMCellType::llmcelltype)
    markers_args <- llmcelltype_args[setdiff(names(llmcelltype_args), llmcelltype_argnames)]
    llmcelltype_args <- llmcelltype_args[intersect(names(llmcelltype_args), llmcelltype_argnames)]

    # Run FindAllMarkers
    log$info("Find the markers for {ident} ...")
    markers_args$object <- sobj
    markers_args$group_by <- ident
    markers_args$ident_1 <- NULL
    markers_args$ident_2 <- NULL
    markers_args$log_prefix <- " * "
    markers_args$log <- log

    sigmarkers <- markers_args$sigmarkers
    markers_args$sigmarkers <- NULL

    llmcelltype_args$input <- do_call(RunSeuratDEAnalysis, markers_args)
    colnames(llmcelltype_args$input)[ncol(llmcelltype_args$input)] <- "cluster"
    if (!is.null(sigmarkers)) {
        log$info("Filtering markers with sigmarkers: {sigmarkers}")
        llmcelltype_args$input <- filter(llmcelltype_args$input, !!parse_expr(sigmarkers))
    }
    rm(sobj)
    rm(markers_args)
    gc()

    # Run LLMCelltype
    log$info("Running LLMCelltype with model '{llmcelltype_args$model}' ...")
    res <- do_call(LLMCellType::llmcelltype, llmcelltype_args)
    if (isTRUE(llmcelltype_args$return_prompt)) {
        stop("LLMCellType prompt:\n\n", res, "\n\n")
    }
    if (is.character(res) && length(res) == 1 && startsWith(res, "Identify cell types")) {
        stop("LLMCellType failed, do you have the correct API key set?")
    }

    # Build mapping (res is named vector: cluster → cell type)
    mapping <- as.list(res)
    list(mapping = mapping)
}
