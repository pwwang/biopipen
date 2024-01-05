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
library(slugify)

setEnrichrSite("Enrichr")

srtfile <- {{ in.srtobj | quote }}
outdir <- {{ out.outdir | quote }}
joboutdir <- {{ job.outdir | quote }}
ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
idents <- {{ envs.idents | r }}
group_by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
subset <- {{ envs.subset | r }}
prefix_each <- {{ envs.prefix_each | r }}
p_adjust <- {{ envs.p_adjust | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
method <- {{ envs.method | r }}
cases <- {{ envs.cases | r: todot = "-" }}

set.seed(8525)

log_info("- Reading Seurat object ...")
srtobj <- readRDS(srtfile)

log_info("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

log_info("- Expanding cases ...")
if (is.null(cases) || length(cases) == 0) {
    cases <- list(
        DEFAULT = list(
            idents = idents,
            group_by = group_by,
            each = each,
            prefix_each = prefix_each,
            p_adjust = p_adjust,
            subset = subset,
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
            subset = subset,
            dbs = dbs,
            sigmarkers = sigmarkers,
            method = method
        )
        cases[[name]] <- case
    }
}

newcases <- list()
sections <- c()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$each)) {
        sections <- c(sections, case$section)
        newcases[[paste0(case$section, ":", name)]] <- case
    } else {
        if (is.null(case$subset)) {
            eachs <- srtobj@meta.data %>% pull(case$each) %>% unique() %>% na.omit()
        } else {
            eachs <- srtobj@meta.data %>% filter(!!parse_expr(case$subset)) %>% pull(case$each) %>% unique() %>% na.omit()
        }
        sections <- c(sections, case$each)
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

# Do enrichment analysis for a case using Enrichr
# Args:
#   case: case name
#   markers: markers dataframe
#   sig: The expression to filter significant markers
do_enrich <- function(info, markers, sig) {
    log_info("  Running enrichment for case: {info$casename}")
    if (nrow(markers) == 0) {
        msg <- paste0("No markers found for case: ", info$casename)
        return(msg)
    }
    markers_sig <- markers %>% filter(!!parse_expr(sig))
    if (nrow(markers_sig) == 0) {
        msg <- paste0("No significant markers found for case: ", info$casename)
        return(msg)
    }
    write.table(
        markers_sig,
        file.path(info$casedir, "markers.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    if (nrow(markers_sig) < 5) {
        msg <- paste0("Too few significant markers found for case: ", info$casename)
        return(msg)
    }

    enriched <- enrichr(markers_sig$gene, dbs)
    for (db in dbs) {
        write.table(
            enriched[[db]],
            file.path(info$casedir, paste0("Enrichr-", db, ".txt")),
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
            file.path(info$casedir, paste0("Enrichr-", db, ".png")),
            res = 100, height = 600, width = 800
        )
        print(
            plotEnrich(enriched[[db]], showTerms = 20, title = db) +
            theme_prism()
        )
        dev.off()
    }
}


do_case <- function(casename) {
    log_info("- Dealing with case: {casename} ...")
    info <- casename_info(casename, create = TRUE)
    case <- cases[[casename]]

    if (sum(!is.na(srtobj@meta.data[[case$group_by]])) == 0) {
        msg = "Not enough cells to run tests."
    } else {
        sobj <- srtobj %>% filter(!is.na(!!sym(case$group_by)))
        if (!is.null(case$subset)) {
            sobj <- srtobj %>% filter(!is.na(!!sym(case$group_by)), !!parse_expr(case$subset))
        }
        df <- tryCatch({
                GetAssayData(sobj, layer = "data")
            }, error = function(e) {
                log_warn("  Error when fetching assay data: {e}")
                NULL
            })
        if (is.null(df)) {
            msg <- "No markers found. May be due to too few cells or features."
        } else {
            genes <- rownames(df)
            # rows: cells, cols: genes
            df <- cbind(as.data.frame(scale(Matrix::t(df))), sobj@meta.data[, case$group_by])
            colnames(df)[ncol(df)] <- "GROUP"

            log_info("  Running tests for case...")
            warn_count <- 0
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
                }, error = function(e) {
                    warn_count <<- warn_count + 1
                    if (warn_count < 10) {
                        log_warn("  Error when testing gene: {gene}")
                        log_warn("  {e}")
                    } else if (warn_count == 10) {
                        log_warn("  Too many errors, will not print more.")
                    }
                    NULL
                })
                if (is.null(res)) {
                    return(NULL)
                }
                res$gene <- gene
                res$method <- case$method
                rownames(res) <- NULL
                res
            }, mc.cores = ncores)
            markers <- do_call(rbind, test_result)
            if (is.null(markers)) {
                msg <- "No markers found. May be due to too few cells."
            } else {
                markers$p_adjust <- p.adjust(markers$p.value, method = case$p_adjust)
                markers <- markers %>% arrange(p_adjust)

                msg <- do_enrich(info, markers, case$sigmarkers)
            }
        }
    }
    if (is.null(msg)) {
        log_info("  Plotting top 10 genes ...")
        markers <- markers %>% head(10)
        plotdir <- file.path(info$casedir, "expr_plots")
        dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

        # Plot the top 10 genes in each group with violin plots
        geneplots <- list()
        for (gene in markers$gene) {
            outfile <- file.path(plotdir, paste0(slugify(gene, tolower = FALSE), ".png"))
            p <- ggplot(df, aes_string(x="GROUP", y=bQuote(gene), fill="GROUP")) +
                geom_violin(alpha = .8) +
                geom_boxplot(width=0.1, fill="white") +
                theme_prism() +
                ylab(paste0("Expression of ", gene))
            png(outfile, res = 100, height = 600, width = 800)
            print(p)
            dev.off()

            geneplots[[length(geneplots) + 1]] <- list(
                kind = "table_image",
                src = outfile,
                name = gene
            )
        }

        add_report(
            list(
                kind = "descr",
                content = paste0(
                    "Top 100 genes selected by ",
                    "<code>", case$method, "</code> across ",
                    "<code>", case$group_by, "</code> and filtered by ",
                    "<code>", html_escape(case$sigmarkers), "</code>"
                )
            ),
            h1 = ifelse(
                info$section == "DEFAULT",
                info$case,
                ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
            ),
            h2 = ifelse(single_section, "Meta-Markers", info$case),
            h3 = ifelse(single_section, "#", "Meta-Markers")
        )
        add_report(
            list(
                name = "Meta-Markers",
                contents = list(list(
                    kind = "table",
                    src = file.path(info$casedir, "markers.txt"),
                    data = list(nrows = 100)
                ))
            ),
            list(
                name = "Volin Plots (Top 10)",
                ui = "table_of_images:4",
                contents = geneplots
            ),
            h1 = ifelse(
                info$section == "DEFAULT",
                info$case,
                ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
            ),
            h2 = ifelse(single_section, "Meta-Markers", info$case),
            h3 = ifelse(single_section, "#", "Meta-Markers"),
            ui = "tabs"
        )
        add_report(
            list(kind = "enrichr", dir = info$casedir),
            h1 = ifelse(
                info$section == "DEFAULT",
                info$case,
                ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
            ),
            h2 = ifelse(single_section, "Enrichment Analysis", info$case),
            h3 = ifelse(single_section, "#", "Enrichment Analysis")
        )
    } else {
        log_warn("  {msg}")
        add_report(
            list(kind = "error", content = msg),
            h1 = ifelse(
                info$section == "DEFAULT",
                info$case,
                ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
            ),
            h2 = ifelse(single_section, "#", info$case)
        )
    }
}

sapply(sort(names(cases)), do_case)
save_report(joboutdir)
