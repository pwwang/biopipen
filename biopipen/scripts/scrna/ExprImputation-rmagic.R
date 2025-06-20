tryCatch(
    {
        # in order to load Rmagic
        workdir <- {{ job.outdir | r }}
        conda_prefix <- Sys.getenv("CONDA_PREFIX")
        setwd(workdir)
        if (!dir.exists("miniconda3")) {
            file.symlink(conda_prefix, "miniconda3")
        }
    },
    error = function(e) {}
)

python <- {{ envs.rmagic_args.python | r }}
Sys.setenv(RETICULATE_PYTHON = ifelse(grepl("/", python, fixed = TRUE), python, Sys.which(python)))
# reticulate::use_python(python, require = TRUE)

library(Rmagic)
Rmagic:::load_pymagic()
pymagic <- tryCatch({
    Rmagic:::pymagic
}, error = function(e) {
    NULL
})
if (is.null(pymagic)) {
    stop("Failed to load pymagic module. Please check your Python environment.\n ",
         "Current python used by reticulate: ", reticulate::py_config()$python)
}

library(Matrix)
library(Seurat)
library(biopipen.utils)

log <- get_logger()

infile <- {{ in.infile | r }}
outfile <- {{ out.outfile | r }}
threshold <- {{ envs.rmagic_args.threshold | r }}

log$info("Loading Seurat object ...")
sobj <- read_obj(infile)

if (threshold > 0) {
    # only use the genes with expression in number of cells greater than threshold
    log$info("Fetching genes with expression great than threshold ({threshold}) ...")
    # get the expression matrix
    layers <- Layers(sobj)
    layer <- ifelse(!"counts" %in% layers, "data", "counts")
    counts <- GetAssayData(sobj, layer = layer)
    # Percent of cells expressing each gene
    dropout_rates <- Matrix::rowSums(counts == 0) / ncol(counts)

    # Genes to impute
    genes_to_impute <- names(dropout_rates[dropout_rates > threshold])

    log$info("- Will impute for {length(genes_to_impute)}/{length(dropout_rates)} genes ...")
    rm(counts)
    rm(dropout_rates)
    gc()
} else {
    genes_to_impute <- NULL
}

# get the expression matrix
data_impute <- t(GetAssayData(sobj, layer = "data"))

log$info("Running MAGIC ...")
check.int.or.null <- function(x) {
    if (is.numeric(x = x)) {
        x <- as.integer(x = x)
    } else if (!is.null(x = x) && is.na(x = x)) {
        x <- NULL
    }
    x
}

check.double.or.null <- function(x) {
    if (is.numeric(x = x)) {
        x <- as.integer(x = x)
    } else if (!is.null(x = x) && is.na(x = x)) {
        x <- NULL
    }
    x
}

check.int.or.string <- function(x, str) {
    if (is.numeric(x = x)) {
        x <- as.integer(x = x)
    } else if (is.null(x = x) || is.na(x = x)) {
        x <- str
    }
    x
}
# the magic function is defined in the Rmagic package
# it has a bug at line 138 when genes are given as a character vector
# See also https://github.com/KrishnaswamyLab/MAGIC/issues/227
magic_patched <- function(
    data,
    genes = NULL,
    knn = 5,
    knn.max = NULL,
    decay = 1,
    t = 3,
    npca = 100,
    solver = "exact",
    init = NULL,
    t.max = 20,
    knn.dist.method = "euclidean",
    verbose = 1,
    n.jobs = 1,
    seed = NULL,
    # deprecated args
    k = NULL, alpha = NULL,
    ...) {
    # check installation
    # if (!reticulate::py_module_available(module = "magic") ||
    #     !exists("pymagic") || is.null(pymagic)) {
    #     Rmagic:::load_pymagic()
    # }
    # check for deprecated arguments
    if (!is.null(k)) {
        message("Argument k is deprecated. Using knn instead.")
        knn <- k
    }
    if (!is.null(alpha)) {
        message("Argument alpha is deprecated. Using decay instead.")
        decay <- alpha
    }
    # validate parameters
    knn <- as.integer(x = knn)
    t.max <- as.integer(x = t.max)
    n.jobs <- as.integer(x = n.jobs)
    npca <- check.int.or.null(npca)
    knn.max <- check.int.or.null(knn.max)
    seed <- check.int.or.null(seed)
    verbose <- check.int.or.null(verbose)
    decay <- check.double.or.null(decay)
    t <- check.int.or.string(t, "auto")
    if (!methods::is(object = data, "Matrix")) {
        data <- as.matrix(x = data)
    }
    # if (length(genes) <= 1 && (is.null(x = genes) || is.na(x = genes))) {
    #                                                  ^^^^^^^^^^^^^^^^ bug here
    if (length(genes) <= 1 && (is.null(x = genes) || (length(genes) == 1 && is.na(x = genes)))) {
        genes <- NULL
        gene_names <- colnames(x = data)
    } else if (is.numeric(x = genes)) {
        gene_names <- colnames(x = data)[genes]
        genes <- as.integer(x = genes - 1)
    } else if (length(x = genes) == 1 && genes == "all_genes") {
        gene_names <- colnames(x = data)
    } else if (length(x = genes) == 1 && genes == "pca_only") {
        gene_names <- paste0("PC", 1:npca)
    } else {
        # character vector
        if (!all(genes %in% colnames(x = data))) {
            warning(paste0(
                "Genes ",
                genes[!(genes %in% colnames(data))],
                " not found.",
                collapse = ", "
            ))
        }
        genes <- which(x = colnames(x = data) %in% genes)
        gene_names <- colnames(x = data)[genes]
        genes <- as.integer(x = genes - 1)
    }
    # store parameters
    params <- list(
        "data" = data,
        "knn" = knn,
        "knn.max" = knn.max,
        "decay" = decay,
        "t" = t,
        "npca" = npca,
        "solver" = solver,
        "knn.dist.method" = knn.dist.method
    )
    # use pre-initialized values if given
    operator <- NULL
    if (!is.null(x = init)) {
        if (!methods::is(init, "magic")) {
            warning("object passed to init is not a phate object")
        } else {
            operator <- init$operator
            operator$set_params(
                knn = knn,
                knn_max = knn.max,
                decay = decay,
                t = t,
                n_pca = npca,
                solver = solver,
                knn_dist = knn.dist.method,
                n_jobs = n.jobs,
                random_state = seed,
                verbose = verbose,
                ...
            )
        }
    }
    if (is.null(x = operator)) {
        operator <- pymagic$MAGIC(
            knn = knn,
            knn_max = knn.max,
            decay = decay,
            t = t,
            n_pca = npca,
            solver = solver,
            knn_dist = knn.dist.method,
            n_jobs = n.jobs,
            random_state = seed,
            verbose = verbose,
            ...
        )
    }
    result <- operator$fit_transform(
        data,
        genes = genes,
        t_max = t.max
    )
    colnames(x = result) <- gene_names
    rownames(x = result) <- rownames(data)
    result <- as.data.frame(x = result)
    result <- list(
        "result" = result,
        "operator" = operator,
        "params" = params
    )
    class(x = result) <- c("magic", "list")
    return(result)
}

data_impute <- magic_patched(data_impute, genes = genes_to_impute)

if (threshold > 0) {
    data <- t(GetAssayData(sobj, layer = "data"))
    data_impute <- cbind(data[, setdiff(colnames(data), genes_to_impute)], Matrix::as.matrix(data_impute$result))
    rm(data)
    gc()
} else {
    # if threshold is 0, then we need to transpose the data back
    data_impute <- t(Matrix::as.matrix(data_impute$result))
}

log$info("Adding imputed data to Seurat object ...")
# Add imputed data to the Seurat object
sobj <- SetAssayData(
    sobj,
    layer = "data",
    new.data = t(data_impute)
)

sobj@misc$impute_method <- "rmagic"

log$info("Saving Seurat object ...")
save_obj(sobj, outfile)
