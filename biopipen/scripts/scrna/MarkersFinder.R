library(dplyr)
library(tibble)
library(Seurat)
library(enrichR)
# library(future)
library(parallel)
setEnrichrSite("Enrichr")

srtobjfile = {{in.srtobj | quote}}
groupfile = {{in.groupfile | quote}}
outdir = {{out.outdir | quote}}
{% if in.casefile %}
cases = {{in.casefile | read | toml_loads | r}}
{% else %}
cases = {{envs.cases | r}}
{% endif %}
dbs = {{envs.dbs | r}}
ncores = {{envs.ncores | r}}

# options(future.globals.maxSize = 80000 * 1024^2)
# plan(strategy = "multicore", workers = ncores_fut)

if (length(cases) == 0) {
    stop("No `envs.cases` specified.")
}

seurat_obj = readRDS(srtobjfile)
# Causing subsets failed to merge
if ("SCT" %in% names(seurat_obj@assays)) {
    seurat_obj[["SCT"]] = NULL
}

groups = read.table(groupfile, row.names=NULL, header=T, sep="\t", check.names = F)
n_samples = ncol(groups) - 1
samples = colnames(groups)[2:(n_samples+1)]
groups = groups %>% rowwise() %>%
    mutate(
        across(2:(n_samples+1),
        ~ strsplit(.x, ";", fixed=TRUE)
    )) %>%
    as.data.frame()

if (is.character(cases) && cases == "ident") {
    cases = list()
    for (ident in unique(Idents(seurat_obj))) {
        cases[[ident]] = list(IDENT=T)
    }
}
if (is.character(cases) && cases == "ALL") {
    cases = list(ALL=list())
}
if ("ALL" %in% names(cases)) {
    cases$ALL = NULL
    allcases = list()
    for (group in unique(as.character(groups[[1]]))) {
        allcases[[group]] = cases
    }
    cases = allcases
}

do_enrich = function(case, markers) {
    casedir = file.path(outdir, case)
    dir.create(casedir, showWarnings = FALSE)
    markers_sig = markers %>%
        filter(p_val_adj < 0.05) %>%
        rownames_to_column("Gene")
    write.table(
        markers_sig,
        file.path(casedir, "markers.txt"),
        sep="\t",
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE
    )

    enriched = enrichr(markers_sig$Gene, dbs)
    for (db in dbs) {
        write.table(
            enriched[[db]],
            file.path(casedir, paste0("Enrichr-", db, ".txt")),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE
        )
        png(
            file.path(casedir, paste0("Enrichr-", db, ".png")),
            res=100, height=1000, width=1000
        )
        if (nrow(markers_sig) == 0) {
            print(ggplot() + annotate("text", x=1, y=1, label="No significant markers."))
        } else {
            print(plotEnrich(enriched[[db]], showTerms = 20, title=db))
        }
        dev.off()
    }
}

do_ident = function(ident) {
    markers = FindMarkers(object = seurat_obj, ident.1 = ident)
    do_enrich(ident, markers)
}

do_failure = function(case, error) {
    casedir = file.path(outdir, case)
    dir.create(casedir, showWarnings = FALSE)
    write.table(
        data.frame(Error=error),
        file.path(casedir, "markers.txt"),
        sep="\t",
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE
    )
}

do_case = function(case) {
    print(paste("- Dealing with case:", case, "..."))
    casepms = cases[[case]]
    if (isTRUE(casepms$IDENT)) {
        do_ident(case)
    } else {
        ident.1 = casepms$ident.1
        ident.2 = casepms$ident.2
        if (is.null(ident.2) && length(unique(as.character(groups[[1]]))) == 2) {
            ident.2 = groups %>%
                filter(.[[1]] != ident.1) %>%
                pull(1) %>%
                as.character() %>%
                unique()
        }
        if (samples == "ALL") {
            case_cells = groups %>%
                filter(.[[1]] == ident.1) %>%
                pull(ALL) %>%
                unlist() %>%
                unname()
        } else {
            case_cells = groups %>%
                filter(.[[1]] == ident.1) %>%
                select(all_of(samples)) %>%
                summarise(across(everything(), c)) %>%
                mutate(across(everything(), ~ list(paste(cur_column(), unlist(.x), sep="_")))) %>%
                unlist() %>%
                unname()
        }
        case_obj = tryCatch({
            subset(seurat_obj, cells = case_cells)
        }, error = function(e) {
            # Not enough cells
            do_failure(case, e$message)
            NULL
        })
        if (is.null(case_obj)) {
            return(NULL)
        }
        case_obj$group = ident.1

        if (is.null(ident.2)) {
            ctrl_obj = subset(seurat_obj, cells = case_cells, invert = TRUE)
            ctrl_obj$group = paste0(ident.1, "_2")
        } else {
            if (samples == "ALL") {
                ctrl_cells = groups %>%
                    filter(.[[1]] == ident.2) %>%
                    pull(ALL) %>%
                    unlist() %>%
                    unname()
            } else {
                ctrl_cells = groups %>% filter(.[[1]] == ident.2) %>%
                    select(all_of(samples)) %>%
                    summarise(across(everything(), c)) %>%
                    mutate(across(everything(), ~ list(paste(cur_column(), unlist(.x), sep="_")))) %>%
                    unlist() %>%
                    unname()
            }
            ctrl_obj = tryCatch({
                subset(seurat_obj, cells = ctrl_cells)
            }, error = function(e) {
                # Not enough cells
                do_failure(case, e$message)
                NULL
            })
            if (is.null(ctrl_obj)) {
                return(NULL)
            }
            ctrl_obj$group = ident.2
        }

        sobj = merge(case_obj, ctrl_obj)
        Idents(sobj) = "group"
        casepms$object = sobj
        markers = do.call(FindMarkers, casepms)
        do_enrich(case, markers)
    }
}

mclapply(names(cases), do_case, mc.cores = ncores)
