# CellTypeAnnotation-scagenttype.R — pure R function, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R

annotate_scagenttype <- function(
    sobj, ident, scagenttype_args, outdir,
    h5ad_path = NULL, case_id = NULL
) {
    log <- get_logger()

    require_package("scagenttype", python = scagenttype_args$python)

    # Use pre-converted h5ad or convert now
    if (is.null(h5ad_path)) {
        case_suffix <- if (!is.null(case_id)) paste0(".", case_id) else ""
        h5ad_path <- file.path(
            outdir, paste0("scagenttype", case_suffix, ".h5ad")
        )
        log$info("Converting Seurat object to h5ad ...")
        ConvertSeuratToAnnData(
            sobj, outfile = h5ad_path,
            assay = scagenttype_args$assay, log = log
        )
    }

    # Output paths
    case_suffix <- if (!is.null(case_id)) paste0("_", case_id) else ""
    cluster_file <- file.path(
        outdir, paste0("scagenttype", case_suffix, ".clusters.tsv")
    )
    config_file <- file.path(
        outdir, paste0("scagenttype", case_suffix, ".config.json")
    )
    outfile <- file.path(
        outdir, paste0("scagenttype", case_suffix, ".txt")
    )

    # Cluster memberships (barcode -> cluster) for the wrapper: the h5ad carries
    # no cluster column, and scAgentType requires one for its per-cluster markers.
    # `ident` is guaranteed to be resolved for cluster-level tools (see main script).
    cluster_df <- data.frame(
        barcode = rownames(sobj@meta.data),
        cluster = as.character(sobj@meta.data[[ident]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    write.table(
        cluster_df, cluster_file,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
    )

    # Config for the wrapper: `envs.scagenttype` minus the R-side and credential
    # keys; `tissue`/`species` are folded into `tissue_context` when not given
    config <- scagenttype_args
    config$python <- NULL
    config$assay <- NULL
    config$api_key <- NULL
    config$base_url <- NULL
    config$tissue <- NULL
    config$species <- NULL
    # Keep the agent's caches (marker DBs, LLM responses) in the job dir
    # instead of the job's working directory (default), which could be inside
    # the source tree. The wrapper clears the LLM-response cache on retry.
    if (is.null(config$cache_dir)) config$cache_dir <- outdir
    if (is.null(config$tissue_context)) {
        tctx <- c(
            if (!is.null(scagenttype_args$tissue)) scagenttype_args$tissue,
            if (!is.null(scagenttype_args$species)) scagenttype_args$species
        )
        if (length(tctx) > 0) {
            config$tissue_context <- paste(tctx, collapse = " / ")
        }
    }
    config <- config[!sapply(config, is.null)]
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("`jsonlite` R package is required for scAgentType annotation")
    }
    writeLines(
        jsonlite::toJSON(config, auto_unbox = TRUE, null = "null", na = "null"),
        config_file
    )

    # Pass credentials to the wrapper process as environment variables, since
    # scAgentType reads the keys from the environment at call time
    child_env <- character(0)
    api <- scagenttype_args$api %||% "openai"
    if (!is.null(scagenttype_args$api_key)) {
        key_env <- c(
            openai = "OPENAI_API_KEY",
            anthropic = "ANTHROPIC_API_KEY",
            google = "GOOGLE_API_KEY"
        )[[api]]
        if (is.na(key_env)) stop(paste0("Unknown api: ", api))
        child_env <- c(child_env, paste0(key_env, "=", scagenttype_args$api_key))
    }
    if (!is.null(scagenttype_args$base_url)) {
        base_env <- c(
            openai = "OPENAI_BASE_URL",
            anthropic = "ANTHROPIC_BASE_URL"
        )[[api]]
        if (is.na(base_env)) {
            stop(paste0("`base_url` is not supported for api: ", api))
        }
        child_env <- c(child_env, paste0(base_env, "=", scagenttype_args$base_url))
    }

    # Run the wrapper script
    biopipen_dir <- get_biopipen_dir(scagenttype_args$python)
    wrapper <- file.path(
        biopipen_dir, "scripts", "scrna", "scagenttype-wrapper.py"
    )
    log$info("Running scAgentType ...")
    rc <- system2(
        scagenttype_args$python,
        c(
            wrapper,
            "-i", h5ad_path,
            "-c", cluster_file,
            "-o", outfile,
            "--config", config_file
        ),
        env = child_env
    )
    if (rc != 0) {
        stop(
            "Failed to run scAgentType. ",
            "Check the job.stderr file to see the error message."
        )
    }

    # Read the per-cell labels and aggregate to a cluster-level mapping.
    # Seurat meta.data columns carry no cell names, so align by `match()` and
    # subset positionally (name-subsetting would return all NAs).
    results <- read.table(outfile, sep = "\t", header = TRUE, row.names = 1)
    idx <- match(rownames(results), rownames(sobj@meta.data))
    if (anyNA(idx)) {
        stop(
            sum(is.na(idx)), " barcodes from scAgentType results are not found ",
            "in the Seurat object. The barcodes were likely changed during ",
            "the h5ad conversion."
        )
    }
    log$info("Aggregating scAgentType results by cluster ...")
    mapping <- majority_vote(
        results[["scagenttype_celltype"]],
        as.character(sobj@meta.data[[ident]][idx])
    )
    list(mapping = mapping)
}
