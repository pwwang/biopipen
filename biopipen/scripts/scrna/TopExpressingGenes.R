source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(tibble)
library(enrichR)
library(rlang)
library(dplyr)

setEnrichrSite("Enrichr")

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
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

print("- Loading Seurat object ...")
srtobj <- readRDS(srtfile)

print("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

print("- Expanding cases ...")
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
for (name in names(cases)) {  # nolint
    case <- cases[[name]]
    if (is.null(case$each) && !is.null(case$ident)) {
        newcases[[paste0(case$section, ":", name)]] <- case
    } else if (is.null(case$each)) {
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
                    key <- paste0(each, kname, ":", ident)
                    if (case$prefix_each) {
                        key <- paste0(case$each, "-", key)
                    }
                    newcases[[key]] <- case
                    newcases[[key]]$ident <- ident
                    newcases[[key]]$group.by <- by  # nolint
                }
            } else {
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

do_enrich <- function(expr, odir) {
    print("  Saving expressions ...")
    write.table(
        expr %>% as.data.frame() %>% rownames_to_column("Gene"),
        file.path(odir, "expr.txt"),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )
    write.table(
        expr %>% as.data.frame() %>% rownames_to_column("Gene") %>% head(n),
        file.path(odir, "exprn.txt"),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE,
        quote = FALSE
    )

    print("  Running enrichment ...")
    enriched <- enrichr(rownames(head(expr, n)), dbs)  # nolint
    for (db in dbs) {
        write.table(
            enriched[[db]],
            file.path(odir, paste0("Enrichr-", db, ".txt")),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
        png(
            file.path(odir, paste0("Enrichr-", db, ".png")),
            res = 100, height = 1000, width = 1000
        )
        print(plotEnrich(enriched[[db]], showTerms = 20, title = db))  # nolint
        dev.off()
    }
}

do_case <- function(casename) {
    print(paste("- Running for case:", casename))
    case <- cases[[casename]]
    parts <- unlist(strsplit(casename, ":"))
    section <- parts[1]
    casename <- paste(parts[-1], collapse = ":")

    print("  Calculating average expression ...")
    avgexpr <- AverageExpression(
        srtobj,
        group.by = case$group.by
    )$RNA[, case$ident, drop = FALSE]
    avgexpr <- avgexpr[order(-avgexpr), , drop = FALSE]

    odir <- file.path(outdir, section, casename)
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    do_enrich(avgexpr, odir)
}

sapply(sort(names(cases)), do_case)
