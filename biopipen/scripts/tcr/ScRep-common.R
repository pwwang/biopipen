suppressPackageStartupMessages({
    library(rlang)
    library(dplyr)
})
#' Common utility functions for single-cell TCR/BCF data analysis/visualization
#'
#'
#' Mutate the data in a Seurat object or a list of data.frames
#' @param screp A Seurat object or a list of data.frames
#' @param mutaters A list of expressions (in characters) to mutate the data
#' @return A Seurat object or a list of data.frames
screp_mutate <- function(screp, mutaters) {
    if (length(mutaters) == 0 || is.null(mutaters)) {
        return (screp)
    }

    if (inherits(screp, "Seurat")) {
        mutaters <- lapply(mutaters, parse_expr)
        screp@meta.data <- mutate(screp@meta.data, !!!mutaters)
    } else {
        screp <- lapply(screp, function(x) {
            mutate(x, !!!mutaters)
        })
    }

    screp
}


#' Filter/subset the data in a Seurat object or a list of data.frames
#'
#' @param screp A Seurat object or a list of data.frames
#' @param subset A character expression to filter the data
#' @return A Seurat object or a list of data.frames
screp_subset <- function(screp, subset) {
    if (inherits(screp, "Seurat")) {
        screp <- eval(parse(text = paste('subset(screp, subset = "', subset, '")')))
    } else {
        screp <- lapply(screp, function(x) {
            dplyr::filter(x, !!parse_expr(subset))
        })
        screp <- screp[sapply(screp, nrow) > 0]
    }

    screp
}
