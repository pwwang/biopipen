suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(immunarch))

#' Expand a Immunarch object into cell-level
#'
#' Here is how the data is expanded:
#' 1. Expand `$data` by Barcode (other columns are copied)
#' 2. Add sample to `Sample` column
#' 3. Add columns from `$meta`
#'
#' @param immdata Immunarch object
#' @return A data frame
#'
#' @example
#' immunarch::immdata$data$MS1 |> head()
#' #   Clones Proportion CDR3.nt            CDR3.aa V.name D.name J.name V.end D.start D.end J.start VJ.ins VD.ins DJ.ins Sequence
#' #    <dbl>      <dbl> <chr>              <chr>   <chr>  <chr>  <chr>  <int>   <int> <int>   <int>  <dbl>  <dbl>  <dbl> <lgl>
#' # 1    539    0.0634  TGTGCCAGCAGCTTACA… CASSLQ… TRBV7… TRBD2  TRBJ2…    14      18    26      29     -1      3      2 NA
#' # 2    320    0.0376  TGTGCCAGCAGCGTGTA… CASSVY… TRBV9  TRBD1  TRBJ2…    13      20    22      29     -1      6      6 NA
#' immunarch::immdata$meta |> head()
#' #   Sample  ID    Sex     Age Status Lane
#' #   <chr>   <chr> <chr> <dbl> <chr>  <chr>
#' # 1 A2-i129 C1    M        11 C      A
#' # 2 A2-i131 C2    M         9 C      A
#' # 3 A2-i133 C4    M        16 C      A
#' # 4 A2-i132 C3    F         6 C      A
#' # 5 A4-i191 C8    F        22 C      B
#' # 6 A4-i192 C9    F        24 C      B
#'
#' @export
expand_immdata <- function(immdata, cell_id = "Barcode") {
    if (!cell_id %in% colnames(immdata$data[[1]])) {
        stop(paste0("cell_id '", cell_id, "' not found in data"))
    }
    do.call(rbind, lapply(names(immdata$data), function(name) {
        # Split barcodes
        dat <- immdata$data[[name]] %>% separate_rows(!!sym(cell_id), sep = ";")
        dat$Sample <- name
        dat <- dat %>% left_join(immdata$meta, by = "Sample", suffix = c("_data", ""))

        dat
    }))
}

#' Filter expanded immdata
#'
#' @param exdata Expanded immdata
#' @param filters Filters
#' @return Filtered data
#'
#' @export
filter_expanded_immdata <- function(exdata, filters, update_clones = FALSE) {
    if (length(filters) == 0) {
        return(exdata)
    }
    out <- exdata %>% dplyr::filter(!!parse_expr(filters))
    if (update_clones) {
        out <- out %>%
            group_by(Sample, CDR3.aa) %>%
            mutate(Clones = n()) %>%
            ungroup() %>%
            group_by(Sample) %>%
            mutate(Proportion = Clones / n()) %>%
            ungroup() %>%
            arrange(Sample, desc(Clones))
    }
    out
}

#' Convert expanded immdata to Immunarch object
#'
#' @param exdata Expanded immdata
#' @param metacols Columns to be added to `$meta`
#' @return Immunarch object
#'
#' @export
immdata_from_expanded <- function(
    exdata,
    metacols = NULL,
    cell_id = "Barcode",
    update_clones = TRUE
) {
    if (is.null(metacols)) {
        metacols = setdiff(colnames(exdata), c(
            "Clones", "Proportion", "CDR3.nt", "CDR3.aa", "V.name", "D.name", "J.name",
            "V.end", "D.start", "D.end", "J.start", "VJ.ins", "VD.ins", "DJ.ins",
            "Sequence", "chain", "Barcode", "raw_clonotype_id", "ContigID", "C.name",
            "CDR1.nt", "CDR2.nt", "CDR1.aa", "CDR2.aa", "FR1.nt", "FR2.nt", "FR3.nt",
            "FR4.nt", "FR1.aa", "FR2.aa", "FR3.aa", "FR4.aa"
        ))
    }
    out <- list(meta = exdata[, metacols, drop = FALSE])
    out$meta <- out$meta[!duplicated(out$meta$Sample), , drop = FALSE]
    out$data <- lapply(
        split(
            exdata[, setdiff(colnames(exdata), metacols), drop = FALSE],
            exdata$Sample
        ),
        function(dat) {
            ncells <- nrow(dat)
            dat_cols <- setdiff(colnames(dat), c("Clones", "Proportion", cell_id))
            dat %>% group_by(CDR3.aa) %>%
                summarise(
                    Clones = ifelse(update_clones, n(), first(Clones)),
                    Proportion = ifelse(update_clones, n() / ncells, first(Proportion)),
                    !!sym(cell_id) := paste0(!!sym(cell_id), collapse = ";"),
                    !!!parse_exprs(sapply(dat_cols, function(x) paste0('first(`', x, '`)'))),
                    .groups = "drop"
                ) %>%
                arrange(desc(Clones))
        }
    )
    out
}
