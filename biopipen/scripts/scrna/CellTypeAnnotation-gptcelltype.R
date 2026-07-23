# CellTypeAnnotation-gptcelltype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_gptcelltype <- function(sobj, ident, gptcelltype_args) {
    library(GPTCelltype)

    log <- get_logger()

    # Set API key
    api_key <- gptcelltype_args$api_key
    gptcelltype_args$api_key <- NULL
    if (is.null(api_key) || api_key == "") {
        stop("`api_key` is required in `gptcelltype_args`")
    }
    Sys.setenv(OPENAI_API_KEY = api_key)

    # Set custom base URL for OpenAI-compatible providers
    base_url <- gptcelltype_args$base_url
    gptcelltype_args$base_url <- NULL
    if (!is.null(base_url) && base_url != "") {
        Sys.setenv(OPENAI_BASE_URL = base_url)
        log$info("Using custom API base URL: {base_url}")
    }

    # Extract model (required), tissuename, and FindAllMarkers args
    model <- gptcelltype_args$model
    gptcelltype_args$model <- NULL
    if (is.null(model) || model == "") {
        stop("`model` is required in `gptcelltype_args`")
    }
    tissuename <- gptcelltype_args$tissuename
    gptcelltype_args$tissuename <- NULL
    assay <- gptcelltype_args$assay
    gptcelltype_args$assay <- NULL

    # Remaining args go to FindAllMarkers
    markers_args <- gptcelltype_args

    # Set identity
    Idents(sobj) <- ident

    # Run FindAllMarkers
    log$info("Running FindAllMarkers ...")
    markers_args$object <- sobj
    markers <- do_call(FindAllMarkers, markers_args)

    # Run GPTCelltype
    log$info("Running GPTCelltype with model '{model}' ...")
    res <- gptcelltype(markers, model = model, tissuename = tissuename)

    # Build mapping (res is named vector: cluster → cell type)
    mapping <- as.list(res)
    list(mapping = mapping)
}
