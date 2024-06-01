suppressPackageStartupMessages({
    library(rlang)
    library(dplyr)
    library(mygene)
})


#@' Convert gene names between different formats
#@'
#@' @param genes A character/integer vector of gene names/ids
#@' @param species A character vector of species names
#@' @param infmt A character vector of input gene name formats
#@'   See the available scopes at
#@'   https://docs.mygene.info/en/latest/doc/data.html#available-fields
#@'   You can use ensg as a shortcut for ensembl.gene
#@' @param outfmt A character vector of output gene name formats
#@' @param dup How to deal with duplicate gene names found.
#@'   "first": keep the first one (default), sorted by score descendingly
#@'   "last": keep the last one, sorted by score descendingly
#@'   "all": keep all of them, each will be a separate row
#@'   "<X>": combine them into a single string, separated by X
#@' @param notfound How to deal with gene names that are not found
#@'   "error": stop with an error message
#@'   "use-query": use the query gene name as the converted gene name
#@'   "skip": skip the gene names that are not found
#@'   "ignore": Same as "skip"
#@'   "na": use NA as the converted gene name (default)
#@' @param suppress_messages Whether to suppress the warning messages
#@' @return A tibble with the query gene names and the converted gene names
#@'   When a gene name is not found, the converted name will be NA
#@'   When duplicate gene names are found, the one with the highest score will be kept
#@' @export
gene_name_conversion <- function(
    genes,
    infmt,
    outfmt,
    dup = "first",
    species = "human",
    notfound = "na",
    suppress_messages = FALSE
) {
    notfound <- arg_match(notfound, c("error", "use-query", "skip", "ignore", "na"))

    if (infmt %in% c("ensg", "ensmusg")) { infmt = "ensembl.gene" }
    if (outfmt %in% c("ensg", "ensmusg")) { outfmt = "ensembl.gene" }

    orig_genes <- genes
    if (infmt == "ensembl.gene") {
        # Remove version numbers from ensembl gene ids
        genes <- gsub("\\..*", "", genes)
    }
    query_df <- tibble(query = genes, orig = orig_genes)

    if (suppress_messages) {
        capture.output(suppressWarnings(suppressMessages({
            out <- queryMany(genes, scopes=infmt, fields=outfmt, species=species) %>%
                as_tibble()
        })))
    } else {
        out <- queryMany(genes, scopes=infmt, fields=outfmt, species=species) %>%
            as_tibble()
    }

    if (nrow(out) == 0) {
        return(tibble(query = orig_genes, converted = NA_character_))
    }

    if (dup == "first") {
        out = out %>% group_by(query) %>% arrange(desc(X_score)) %>%
            slice_head(n=1) %>% ungroup() %>% dplyr::select(all_of(c("query", outfmt)))
    } else if (dup == "last") {
        out = out %>% group_by(query) %>% arrange(X_score) %>%
            slice_head(n=1) %>% ungroup() %>% dplyr::select(all_of(c("query", outfmt)))
    } else if (dup != "all") {
        out = out %>% group_by(query) %>% arrange(desc(X_score)) %>%
            summarise(!!sym(outfmt) := paste(unique(!!sym(outfmt)), collapse=dup))
    }
    out <- query_df %>%
        left_join(out, by="query") %>%
        dplyr::select(-"query") %>%
        dplyr::select(query = orig, everything())

    if (notfound == "error") {
        if (any(is.na(out[[outfmt]]))) {
            nagenes = out %>% filter(is.na(.[[outfmt]])) %>% pull("query")
            stop(paste("Query genes not found:", paste(nagenes, collapse=",")))
        }
    } else if (notfound == "use-query") {
        out = out %>% mutate(!!sym(outfmt) := coalesce(!!sym(outfmt), query))
    } else if (notfound == "skip" || notfound == "ignore") {
        out = out %>% filter(!is.na(!!sym(outfmt)))
    }

    return(out)
}
