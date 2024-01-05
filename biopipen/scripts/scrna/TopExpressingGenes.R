source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(tibble)
library(enrichR)
library(rlang)
library(dplyr)
library(slugify)
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
cases <- {{ envs.cases | r: todot = "-" }}  # nolint

set.seed(8525)

log_info("Loading Seurat object ...")
srtobj <- readRDS(srtfile)

log_info("Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

log_info("Expanding cases ...")
if (is.null(cases) || length(cases) == 0) {
    cases <- list(
        DEFAULT = list(
            ident = ident,
            group.by = group.by,
            each = each,
            prefix_each = prefix_each,
            section = section,
            dbs = dbs,
            n = n
        )
    )
} else {
    cases <- lapply(cases, function(cs) {
        list_setdefault(
            cs,
            ident = ident,
            group.by = group.by,
            each = each,
            prefix_each = prefix_each,
            section = section,
            dbs = dbs,
            n = n
        )
    })
}

# Expand each and ident
newcases <- list()
sections <- c()
for (name in names(cases)) {  # nolint
    case <- cases[[name]]
    if (is.null(case$each) && !is.null(case$ident)) {
        sections <- c(sections, case$section)
        newcases[[paste0(case$section, ":", name)]] <- case
    } else if (is.null(case$each)) {
        sections <- c(sections, name)
        idents <- srtobj@meta.data %>%
            pull(case$group.by) %>%
            unique() %>%
            na.omit()
        for (ident in idents) {
            key <- paste0(name, ":", ident)
            newcases[[key]] <- case
            newcases[[key]]$ident <- ident
        }
    } else {
        eachs <- srtobj@meta.data %>% pull(case$each) %>% unique() %>% na.omit()
        for (each in eachs) {
            by <- make.names(paste0(".", name, "_", each))
            srtobj@meta.data <- srtobj@meta.data %>% mutate(
                !!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group.by),
                    NA
                )
            )
            if (is.null(case$ident)) {
                idents <- srtobj@meta.data %>%
                    pull(case$group.by) %>%
                    unique() %>%
                    na.omit()
                for (ident in idents) {
                    kname <- if (name == "DEFAULT") "" else paste0("-", name)
                    sections <- c(sections, paste0(each, kname))
                    key <- paste0(each, kname, ":", ident)
                    if (case$prefix_each) {
                        key <- paste0(
                            ifelse(case$each == "seurat_clusters", "Cluster", case$each),
                            " - ",
                            key
                        )
                    }
                    newcases[[key]] <- case
                    newcases[[key]]$ident <- ident
                    newcases[[key]]$group.by <- by  # nolint
                }
            } else {
                sections <- c(sections, case$each)
                key <- paste0(case$each, ":", each)
                if (name != "DEFAULT") {
                    key <- paste0(key, " - ", name)
                }
                newcases[[key]] <- case
            }
        }
    }
}
cases <- newcases
single_section <- length(unique(sections)) == 1

casename_info <- function(casename, create = FALSE) {
    sec_case_names <- strsplit(casename, ":")[[1]]
    cname <- paste(sec_case_names[-1], collapse = ":")

    out <- list(
        casename = casename,
        section = sec_case_names[1],
        case = cname,
        section_slug = slugify(sec_case_names[1], tolower = FALSE),
        case_slug = slugify(cname, tolower = FALSE)
    )
    out$casedir <- file.path(outdir, out$section_slug, out$case_slug)
    if (create) {
        dir.create(out$casedir, showWarnings = FALSE, recursive = TRUE)
    }
    out
}

do_enrich <- function(expr, odir) {
    log_info("  Saving expressions ...")
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

    log_info("  Running enrichment ...")
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
            log_info(paste0("  No enriched terms for ", db))
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
    info <- casename_info(casename, create = TRUE)

    log_info("  Calculating average expression ...")
    assay <- DefaultAssay(srtobj)
    avgexpr <- AverageExpression(
        srtobj,
        group.by = case$group.by,
        assays = assay
    )[[assay]]
    # https://github.com/satijalab/seurat/issues/7893
    colnames(avgexpr) <- as.character(unique(srtobj@meta.data[[case$group.by]]))
    avgexpr <- avgexpr[, case$ident, drop = FALSE]
    avgexpr <- avgexpr[order(-avgexpr), , drop = FALSE]

    do_enrich(avgexpr, info$casedir)

    add_case_report(info)
}

add_case_report <- function(info) {
    log_info("  Adding case report ...")
    h1 = ifelse(
        info$section == "DEFAULT",
        info$case,
        ifelse(
            single_section,
            paste0(
                ifelse(info$section == "seurat_clusters", "Cluster", info$section),
                " - ",
                info$case
            ),
            info$section
        )
    )
    h2 = ifelse(
        info$section == "DEFAULT",
        "#",
        ifelse(single_section, "#", info$case)
    )

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

sapply(sort(names(cases)), do_case)
save_report(joboutdir)
