source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(Seurat)
library(enrichR)
library(ggplot2)
library(ggprism)
library(parallel)
library(tidyseurat)

setEnrichrSite("Enrichr")

srtfile <- {{ in.srtobj | quote }}
outdir <- {{ out.outdir | quote }}
ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
idents <- {{ envs.idents | r }}
group_by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
p_adjust <- {{ envs.p_adjust | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
method <- {{ envs.method | r }}
cases <- {{ envs.cases | r: todot = "-" }}

set.seed(8525)

print("- Reading Seurat object ...")
srtobj <- readRDS(srtfile)

print("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

print("- Expanding cases ...")
if (is.null(cases) || length(cases) == 0) {
    cases <- list(
        DEFAULT = list(
            idents = idents,
            group_by = group_by,
            each = each,
            prefix_each = prefix_each,
            p_adjust = p_adjust,
            section = section,
            dbs = dbs,
            sigmarkers = sigmarkers,
            method = method
        )
    )
} else {
    for (name in names(cases)) {
        case <- list_setdefault(
            cases[[name]],
            idents = idents,
            group_by = group_by,
            each = each,
            prefix_each = prefix_each,
            p_adjust = p_adjust,
            section = section,
            dbs = dbs,
            sigmarkers = sigmarkers,
            method = method
        )
        cases[[name]] <- case
    }
}

newcases <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$each)) {
        newcases[[paste0(case$section, ":", name)]] <- case
    } else {
        eachs <- srtobj@meta.data %>% pull(case$each) %>% unique() %>% na.omit()
        for (each in eachs) {
            by = make.names(paste0(".", name, "_", case$each, "_", each))
            idents <- case$idents
            if (is.null(idents) || length(idents) == 0) {
                srtobj@meta.data = srtobj@meta.data %>%
                    mutate(
                        !!sym(by) := if_else(!!sym(case$each) == each, !!sym(case$group_by), NA)
                    )
                idents <- srtobj@meta.data %>% pull(case$group_by) %>% unique() %>% na.omit()
            } else {
                srtobj@meta.data = srtobj@meta.data %>%
                    mutate(
                        !!sym(by) := if_else(
                            !!sym(case$each) == each & !!sym(case$group_by) %in% case$idents,
                            !!sym(case$group_by),
                            NA
                        )
                    )
            }

            key <- paste0(case$each, ":", each)
            if (name != "DEFAULT") {
                key <- paste0(key, " - ", name)
            }
            newcases[[key]] <- case
            newcases[[key]]$group_by <- by
            newcases[[key]]$idents <- idents
        }
    }
}
cases <- newcases


# Do enrichment analysis for a case using Enrichr
# Args:
#   case: case name
#   markers: markers dataframe
#   sig: The expression to filter significant markers
do_enrich <- function(case, markers, sig) {
    print(paste("  Running enrichment for case:", case))
    parts <- strsplit(case, ":")[[1]]
    sec <- parts[1]
    case <- paste0(parts[-1], collapse = ":")
    casedir <- file.path(outdir, sec, case)
    dir.create(casedir, showWarnings = FALSE, recursive = TRUE)
    if (nrow(markers) == 0) {
        print(paste("  No markers found for case:", case))
        cat("No markers found.", file = file.path(casedir, "error.txt"))
        return()
    }
    markers_sig <- markers %>% filter(!!parse_expr(sig))
    if (nrow(markers_sig) == 0) {
        print(paste("  No significant markers found for case:", case))
        cat("No significant markers.", file = file.path(casedir, "error.txt"))
        return()
    }
    write.table(
        markers_sig,
        file.path(casedir, "markers.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    if (nrow(markers_sig) < 5) {
        for (db in dbs) {
            write.table(
                data.frame(Warning = "Not enough significant markers."),
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res = 100, height = 200, width = 1000
            )
            print(
                ggplot() +
                    annotate(
                        "text",
                        x = 1,
                        y = 1,
                        label = "Not enough significant markers."
                    ) +
                    theme_classic()
            )
            dev.off()
        }
    } else {
        enriched <- enrichr(markers_sig$gene, dbs)
        for (db in dbs) {
            write.table(
                enriched[[db]],
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res = 100, height = 1000, width = 1000
            )
            print(plotEnrich(enriched[[db]], showTerms = 20, title = db))
            dev.off()
        }
    }
}


do_case <- function(casename) {
    cat(paste("- Dealing with case:", casename, "...\n"))
    case <- cases[[casename]]
    sobj <- srtobj %>% tidyseurat::filter(!is.na(!!sym(case$group_by)))
    df <- GetAssayData(sobj, slot = "data", assay = "RNA")
    genes <- rownames(df)
    # rows: cells, cols: genes
    df <- cbind(as.data.frame(scale(Matrix::t(df))), sobj@meta.data[, case$group_by])
    colnames(df)[ncol(df)] <- "GROUP"

    cat(paste("  Running tests for case...\n"))
    test_result <- mclapply(genes, function(gene) {
        fm <- as.formula(paste(bQuote(gene), "~ GROUP"))
        res <- tryCatch({
            if (case$method == "anova") {
                r <- summary(aov(fm, data = df))[[1]]
                data.frame(
                    statistic = r[1, "F value"],
                    p.value = r[1, "Pr(>F)"],
                    sumsq = r[1, "Sum Sq"],
                    meansq = r[1, "Mean Sq"]
                )
            } else {
                r <- kruskal.test(fm, data = df)
                data.frame(statistic = r$statistic, p.value = r$p.value)
            }
        }, error = function(e) NULL)
        if (is.null(res)) {
            return(NULL)
        }
        res$gene <- gene
        res$method <- case$method
        rownames(res) <- NULL
        res
    }, mc.cores = ncores)
    markers <- do_call(rbind, test_result)
    markers$p_adjust <- p.adjust(markers$p.value, method = case$p_adjust)
    markers <- markers %>% arrange(p_adjust)
    do_enrich(casename, markers, case$sigmarkers)

    print(paste("  Plotting top 10 genes ...\n"))
    markers <- markers %>% head(10)
    parts <- strsplit(casename, ":")[[1]]
    sec <- parts[1]
    casename <- paste0(parts[-1], collapse = ":")
    plotdir <- file.path(outdir, sec, casename, "plots")
    dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

    # Plot the top 10 genes in each group with violin plots
    for (gene in markers$gene) {
        outfile = file.path(plotdir, paste0(gene, ".png"))
        p = ggplot(df, aes_string(x="GROUP", y=bQuote(gene), fill="GROUP")) +
            geom_violin(alpha = .8) +
            geom_boxplot(width=0.1, fill="white") +
            theme_prism() +
            ylab(paste0("Expression of ", gene))
        png(outfile, res = 100, height = 800, width = 1000)
        print(p)
        dev.off()
    }
}

sapply(sort(names(cases)), do_case)
