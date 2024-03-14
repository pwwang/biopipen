source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(tibble)
library(enrichR)
library(rlang)
library(dplyr)
library(ggprism)

setEnrichrSite("Enrichr")

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
mutaters <- {{ envs.mutaters | r }}
ident <- {{ envs.ident | r }}
group.by <- {{ envs["group-by"] | r }}  # nolint
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
n <- {{ envs.n | r }}
sset <- {{ envs.subset | r }}
cases <- {{ envs.cases | r: todot = "-" }}  # nolint

set.seed(8525)

log_info("- Loading Seurat object ...")
srtobj <- readRDS(srtfile)
assay <- DefaultAssay(srtobj)

log_info("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    ident = ident,
    group.by = group.by,
    each = each,
    prefix_each = prefix_each,
    section = section,
    dbs = dbs,
    n = n,
    subset = sset
)

expand_each <- function(name, case) {
    outcases <- list()
    no_each <- is.null(case$each) || nchar(case$each) == 0
    no_ident <- is.null(case$ident)
    has_section <- !is.null(case$section) && case$section != "DEFAULT"
    if (no_each && !no_ident) {
        # single cases
        if (is.null(case$section) || case$section == "DEFAULT") {
            outcases[[name]] <- case
        } else {
            outcases[[paste0(case$section, "::", name)]] <- case
        }
    } else if (no_each) {  # no_ident
        # expanding idents
        if (has_section) {
            log_warn("  Ignoring `section` in case `{name}` when no `ident` is set.")
            case$section <- NULL
        }
        if (!is.null(case$subset)) {
            idents <- srtobj@meta.data %>% filter(!!parse_expr(case$subset)) %>%
                pull(case$group.by) %>% unique() %>% na.omit() %>% as.vector()
        } else {
            idents <- srtobj@meta.data %>%
                pull(case$group.by) %>% unique() %>% na.omit() %>% as.vector()
        }

        for (ident in idents) {
            key <- paste0(name, "::", ident)
            outcases[[key]] <- case
            outcases[[key]]$ident <- ident
            outcases[[key]]$section <- name
        }
    } else {  # has_each
        if (no_ident) {
            stop("  `ident` must be set when `each` is set for case `{name}`.")
        }
        # expanding eachs
        if (has_section) {
            log_warn("  Ignoring `section` in case `{name}` when `each` is set.")
            case$section <- NULL
        }

        if (!is.null(case$subset)) {
            eachs <- srtobj@meta.data %>% filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% unique() %>% na.omit() %>% as.vector()
        } else {
            eachs <- srtobj@meta.data %>%
                pull(case$each) %>% unique() %>% na.omit() %>% as.vector()
        }

        for (each in eachs) {
            by <- make.names(paste0(".", name, "_", case$each,"_", each))
            srtobj@meta.data <<- srtobj@meta.data %>% mutate(
                !!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group.by),
                    NA
                )
            )

            if (isTRUE(case$prefix_each)) {
                key <- paste0(name, "::", case$each, " - ", each)
            } else {
                key <- paste0(name, "::", each)
            }
            outcases[[key]] <- case
            outcases[[key]]$section <- name
            outcases[[key]]$group.by <- by
        }
    }
    outcases
}

log_info("- Expanding cases ...")
cases <- expand_cases(cases, defaults, expand_each)

do_enrich <- function(expr, odir) {
    log_debug("  Saving expressions ...")
    expr <- expr %>% as.data.frame()
    colnames(expr) <- c("Expression")
    expr <- expr %>% rownames_to_column("Gene") %>% select(Gene, Expression)
    write.table(
        expr,
        file.path(odir, "expr.txt"),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )
    write.table(
        expr %>% head(n),
        file.path(odir, "exprn.txt"),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )

    log_debug("  Running enrichment ...")
    enriched <- enrichr(head(expr$Gene, n), dbs)  # nolint
    for (db in dbs) {
        write.table(
            enriched[[db]],
            file.path(odir, paste0("Enrichr-", db, ".txt")),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )

        if (nrow(enriched[[db]]) == 0) {
            log_warn(paste0("  No enriched terms for ", db))
            next
        }

        png(
            file.path(odir, paste0("Enrichr-", db, ".png")),
            res = 100, height = 1000, width = 1000
        )
        print(
            plotEnrich(enriched[[db]], showTerms = 20, title = db) +
            theme_prism()
        )
        dev.off()
    }
}

do_case <- function(casename) {
    log_info("- Running for case: {casename} ...")
    case <- cases[[casename]]
    info <- casename_info(casename, cases, outdir, create = TRUE)

    log_debug("  Calculating average expression ...")
    if (!is.null(case$subset)) {
        tryCatch({
            sobj <- subset(srtobj, !!parse_expr(case$subset))
        }, error = function(e) {
            log_warn("  No cells found for the subset, skipping ...")
        })
    } else {
        sobj <- srtobj
    }
    avgexpr <- AverageExpression(
        sobj,
        group.by = case$group.by,
        assays = assay
    )[[assay]]
    # https://github.com/satijalab/seurat/issues/7893
    colnames(avgexpr) <- as.character(unique(sobj@meta.data[[case$group.by]]))
    avgexpr <- avgexpr[, case$ident, drop = FALSE]
    avgexpr <- avgexpr[order(-avgexpr), , drop = FALSE]

    do_enrich(avgexpr, info$casedir)

    add_case_report(info)
}

add_case_report <- function(info) {
    log_debug("  Adding case report ...")
    h1 = info$h1
    h2 = info$h2

    if (!is.null(info$error)) {
        add_report(
            list(
                kind = "descr",
                content = paste0("Top ", n, " expressing genes")
            ),
            list(kind = "error", content = info$error),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Top Expressing Genes", h2),
            h3 = ifelse(h2 == "#", "#", "Top Expressing Genes")
        )
    } else {
        add_report(
            list(
                kind = "descr",
                content = paste0("Top ", n, " expressing genes")
            ),
            list(
                kind = "table",
                src = file.path(info$casedir, "exprn.txt")
            ),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Top Expressing Genes", h2),
            h3 = ifelse(h2 == "#", "#", "Top Expressing Genes")
        )

        add_report(
            list(
                kind = "descr",
                content = paste0("Enrichment analysis for the top ", n, " expressing genes")
            ),
            list(kind = "enrichr", dir = info$casedir),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Enrichment Analysis", h2),
            h3 = ifelse(h2 == "#", "#", "Enrichment Analysis")
        )
    }
}

sapply(sort(names(cases)), do_case)
save_report(joboutdir)
