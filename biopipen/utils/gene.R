library(mygene)
library(dplyr)

gene_name_conversion = function(
    genes,
    species,
    infmt,
    outfmt,
    notfound
) {
    out = queryMany(
        genes,
        scopes=infmt,
        fields=outfmt,
        species=species
    ) %>% as.data.frame() %>% group_by(
        query
    ) %>% arrange(
        desc(X_score)
    ) %>% slice_head(n=1) %>% select(
        -c(X_id, X_score)
    )

    if ("notfound" %in% colnames(out)) {
        out = out %>% select(-c("notfound"))
    }

    if (length(outfmt) == 1 && "," %in% outfmt) {
        outfmt = trimws(unlist(strsplit(outfmt, ",", fixed=TRUE)))
    }

    out = tibble(query=genes) %>% left_join(out, by="query")
    if (notfound == "use-query") {
        out = out %>% mutate(
            across(
                outfmt,
                function(col, query) if_else(is.na(col), query, col),
                query=query
            )
        )
    } else if (notfound == "error" && any(is.na(out[[outfmt[1]]]))) {
        nagenes = out %>% filter(is.na(.[[outfmt[1]]])) %>% pull("query")
        stop(paste("Query genes not found:", paste(nagenes, collapse=",")))
    } else if (notfound == "skip") {
        out = out %>% filter(!is.na(.[[outfmt[1]]]))
    }

    return out
}
