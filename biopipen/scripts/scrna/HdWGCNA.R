library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(tidyseurat)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(patchwork)
library(biopipen.utils)

log <- get_logger()
reporter <- get_reporter()

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
ncores <- {{envs.ncores | r}}
mutaters <- {{envs.mutaters | r}}
subset <- {{envs.subset | r}}
seed <- {{envs.seed | r}}
cache <- {{envs.cache | r}}
use_pseudobulk <- {{envs.use_pseudobulk | r}}
use_consensus <- {{envs.use_consensus | r}}
ref_srtfile <- {{envs.ref_srtobj | r}}

qs2::qopt("nthreads", value = ncores)
set.seed(seed)
options(bitmapType = "cairo")
# hdWGCNA's future.apply calls can export multi-GB globals (e.g. the TOM
# matrix in ModuleUMAPPlot); raise the default 500MiB cap on every plan
options(future.globals.maxSize = Inf)
if (ncores > 1) {
    future::plan(future::multicore, workers = ncores, maxSizeOfObjects = Inf)
    WGCNA::enableWGCNAThreads(nThreads = ncores)
}

if (isTRUE(cache)) { cache <- joboutdir }

log$info("Loading the Seurat object ...")
srtobj <- read_obj(srtfile)

log$info("Applying mutaters if any ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj <- MutateSeuratMeta(srtobj, mutaters, log = log)
}

if (!is.null(subset)) {
    log$info("Subsetting the cells ...")
    srtobj <- srtobj %>% filter(!!parse_expr(subset))
}

ref_srtobj <- NULL
ref_is_self <- FALSE
if (!is.null(ref_srtfile) && !identical(ref_srtfile, "self")) {
    # resolve relative paths against the directory of the input object
    if (!file.exists(ref_srtfile)) {
        ref_srtfile <- file.path(dirname(srtfile), ref_srtfile)
    }
    log$info(glue("Loading the reference Seurat object from {ref_srtfile} ..."))
    ref_srtobj <- read_obj(ref_srtfile)
} else if (identical(ref_srtfile, "self")) {
    # use the query object itself as the reference (self-preservation test)
    ref_is_self <- TRUE
}

############## helper functions ##############
hd_formals <- function(fn) names(formals(get(fn, asNamespace("hdWGCNA"))))

hd_args <- function(args, fn) {
    fmls <- hd_formals(fn)
    if ("..." %in% fmls) args else args[intersect(names(args), fmls)]
}

save_plot_list <- function(plots, prefix, devpars, ncol = 2, formats = c("png")) {
    # hdWGCNA functions like ModuleRadarPlot / PlotKMEs already return a
    # combined patchwork (with combine=TRUE / ncol defaults); re-wrapping
    # it nests the patchwork and squeezes every panel into the first column
    if (inherits(plots, "patchwork")) {
        save_plot(plots, prefix, devpars, formats = formats)
    } else if (length(plots) == 1) {
        save_plot(plots[[1]], prefix, devpars, formats = formats)
    } else {
        save_plot(patchwork::wrap_plots(plots, ncol = ncol), prefix, devpars, formats = formats)
    }
}

# report sections with a single case (e.g. "Soft Powers") get an h2
# identical to the h1; drop the redundant h2 in that case. Also drop it
# when the case name is only a plural variant of the section (e.g.
# "Module Radar" under the "Module Radars" section).
hs_section <- function(section, name = NULL) {
    if (
        is.null(name) ||
        identical(section, name) ||
        identical(paste0(section, "s"), name) ||
        identical(section, paste0(name, "s"))
    ) {
        section
    } else {
        c(section, name)
    }
}

# fill in missing device dimensions from the number of panels, so that
# multi-group/module plots are not squeezed into the default 8x6in box
# NOTE: biopipen.utils::save_plot takes width/height in PIXELS (defaults
# 800x600); inches would create near-empty devices (and patchwork errors
# with "The viewport may be too small to show this patchwork")
auto_devpars <- function(devpars, npanels, ncol = 1, per = 2.5) {
    devpars <- devpars %||% list()
    if (npanels > 1) {
        res <- devpars$res %||% 100
        if (is.null(devpars$width)) {
            devpars$width <- max(800, (min(npanels, ncol) * per + 1) * res)
        }
        if (is.null(devpars$height)) {
            devpars$height <- max(600, (ceiling(npanels / ncol) * per + 1) * res)
        }
    }
    devpars
}

# merge PDFs (one per module) into a single file with ghostscript, so the
# report embeds one PDF per plot case instead of dozens. Returns FALSE when
# gs is unavailable; the caller then falls back to the individual files.
merge_pdfs <- function(pdfs, out) {
    gs <- Sys.which("gs")
    if (!nzchar(gs) || length(pdfs) < 2) return(FALSE)
    code <- suppressWarnings(system2(gs, c(
        "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.4",
        paste0("-sOutputFile=", shQuote(out)),
        shQuote(pdfs)
    )))
    identical(code, 0L)
}

tables_dir <- file.path(outdir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

############## core ##############
SetupForWGCNAArgs <- {{envs.SetupForWGCNA | r: todot="-"}}
MetacellsByGroupsArgs <- {{envs.MetacellsByGroups | r: todot="-"}}
NormalizeMetacellsArgs <- {{envs.NormalizeMetacells | r: todot="-"}}
ScaleMetacellsArgs <- {{envs.ScaleMetacells | r: todot="-"}}
RunPCAMetacellsArgs <- {{envs.RunPCAMetacells | r: todot="-"}}
RunHarmonyMetacellsArgs <- {{envs.RunHarmonyMetacells | r: todot="-"}}
RunUMAPMetacellsArgs <- {{envs.RunUMAPMetacells | r: todot="-"}}
AggregatePseudobulkArgs <- {{envs.AggregatePseudobulk | r: todot="-"}}
NormalizeCountsArgs <- {{envs.NormalizeCounts | r: todot="-"}}
SetDatExprArgs <- {{envs.SetDatExpr | r: todot="-"}}
TestSoftPowersArgs <- {{envs.TestSoftPowers | r: todot="-"}}
SetMultiExprArgs <- {{envs.SetMultiExpr | r: todot="-"}}
TestSoftPowersConsensusArgs <- {{envs.TestSoftPowersConsensus | r: todot="-"}}
ConstructNetworkArgs <- {{envs.ConstructNetwork | r: todot="-"}}
ModuleEigengenesArgs <- {{envs.ModuleEigengenes | r: todot="-"}}
ModuleConnectivityArgs <- {{envs.ModuleConnectivity | r: todot="-"}}
ResetModuleNamesArgs <- {{envs.ResetModuleNames | r: todot="-"}}
ResetModuleColorsArgs <- {{envs.ResetModuleColors | r: todot="-"}}
GetHubGenesArgs <- {{envs.GetHubGenes | r: todot="-"}}
ModuleExprScoreArgs <- {{envs.ModuleExprScore | r: todot="-"}}

wgcna_name <- SetupForWGCNAArgs$wgcna_name %||% "biopipen"
SetupForWGCNAArgs$wgcna_name <- NULL

{% include biopipen_dir + "/scripts/scrna/HdWGCNA-core.R" %}

############## plots ##############
plots_defaults <- {{envs.plots_defaults | r: todot="-"}}
plots <- {{envs.plots | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/HdWGCNA-plots.R" %}

############## downstream ##############
dmes_defaults <- {{envs.dmes_defaults | r: todot="-"}}
dmes <- {{envs.dmes | r: todot="-", skip=1}}
module_trait_corr_defaults <- {{envs.module_trait_corr_defaults | r: todot="-"}}
module_trait_corr <- {{envs.module_trait_corr | r: todot="-", skip=1}}
enrich_defaults <- {{envs.enrich_defaults | r: todot="-"}}
enrich <- {{envs.enrich | r: todot="-", skip=1}}
gsea_defaults <- {{envs.gsea_defaults | r: todot="-"}}
gsea <- {{envs.gsea | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/HdWGCNA-downstream.R" %}

############## preservation ##############
module_preservations_defaults <- {{envs.module_preservations_defaults | r: todot="-"}}
module_preservations <- {{envs.module_preservations | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/HdWGCNA-preservation.R" %}

############## tf network ##############
tf_network <- {{envs.tf_network | r: todot="-"}}
{% include biopipen_dir + "/scripts/scrna/HdWGCNA-tf-network.R" %}

############## motifs ##############
motifs <- {{envs.motifs | r: todot="-"}}
{% include biopipen_dir + "/scripts/scrna/HdWGCNA-motifs.R" %}

############## tables ##############
log$info("Writing tables ...")

mod_table <- tryCatch({
    do_call("GetModules", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
}, error = function(e) {
    log$warn(glue("Failed to get the module table: {conditionMessage(e)}"))
    NULL
})
if (!is.null(mod_table)) {
    write_table(mod_table, file.path(tables_dir, "modules.tsv"))
}

hub_genes <- tryCatch({
    do_call("GetHubGenes", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), GetHubGenesArgs))
}, error = function(e) {
    log$warn(glue("Failed to get the hub genes: {conditionMessage(e)}"))
    NULL
})
if (!is.null(hub_genes)) {
    write_table(hub_genes, file.path(tables_dir, "hub_genes.tsv"))
}

power_table <- tryCatch({
    do_call("GetPowerTable", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
}, error = function(e) {
    log$warn(glue("Failed to get the power table: {conditionMessage(e)}"))
    NULL
})
if (!is.null(power_table)) {
    write_table(power_table, file.path(tables_dir, "power_table.tsv"))
}

MEs <- tryCatch({
    do_call("GetMEs", list(seurat_obj = srtobj, wgcna_name = wgcna_name))
}, error = function(e) {
    log$warn(glue("Failed to get the module eigengenes: {conditionMessage(e)}"))
    NULL
})
if (!is.null(MEs)) {
    write_table(MEs, file.path(tables_dir, "MEs.tsv"), row.names = TRUE)
}

############## report: data ##############
# the hub genes table is reported next to the Module UMAP plot when one is
# configured; only add it here (Data/Tables) as a fallback otherwise
has_module_umap <- any(vapply(plots, function(x) identical(x$kind, "module_umap"), logical(1)))
data_entries <- list(
    list(name = "Modules", contents = list(
        list(kind = "descr", content = "The genes and their module assignments."),
        list(kind = "table", src = file.path(tables_dir, "modules.tsv"), data = list(nrows = 100))
    ))
)
if (!has_module_umap) {
    data_entries <- c(data_entries, list(list(name = "Hub Genes", contents = list(
        list(kind = "descr", content = "The hub genes of the modules."),
        list(kind = "table", src = file.path(tables_dir, "hub_genes.tsv"), data = list(nrows = 100))
    ))))
}
data_entries <- c(data_entries, list(
    list(name = "Power Table", contents = list(
        list(kind = "descr", content = "The soft power selection table."),
        list(kind = "table", src = file.path(tables_dir, "power_table.tsv"), data = list(nrows = 100))
    )),
    list(name = "Module Eigengenes", contents = list(
        list(kind = "descr", content = "The module eigengene values of the cells."),
        list(kind = "table", src = file.path(tables_dir, "MEs.tsv"), data = list(nrows = 100))
    ))
))
do.call(reporter$add2, c(data_entries, list(hs = c("Data", "Tables"), ui = "tabs")))

############## report: introduction ##############
mode_str <- if (use_consensus) {
    "a consensus network"
} else if (use_pseudobulk) {
    "a pseudobulk network"
} else {
    "a metacell network"
}

soft_power_info <- tryCatch({
    params <- GetWGCNAParams(srtobj, wgcna_name)
    if (!is.null(params$power)) {
        glue("The soft power of {params$power} was selected for the network construction.")
    } else {
        "The soft power was automatically selected by the analysis."
    }
}, error = function(e) "The soft power was automatically selected by the analysis.")

# `module` is a factor — compare labels, not level codes (see HdWGCNA-plots.R)
n_modules <- if (!is.null(mod_table)) sum(as.character(mod_table$module) != "grey") else NA
n_genes <- if (!is.null(mod_table)) nrow(mod_table) else NA

metacell_info <- if (use_pseudobulk) {
    glue(
        "The cells were aggregated into pseudobulk samples by `AggregatePseudobulk` ",
        "(replicate column: `{AggregatePseudobulkArgs$replicate_col %||% 'Sample'}`, ",
        "group column: `{AggregatePseudobulkArgs$group_col %||% 'seurat_clusters'}`) ",
        "and normalized with `NormalizeCounts` (method: `{NormalizeCountsArgs$method %||% 'VST'}`)."
    )
} else {
    glue(
        "The metacells were constructed by `MetacellsByGroups` ",
        "(group by: `{paste(MetacellsByGroupsArgs$group.by, collapse = ', ')}`, ",
        "ident group: `{MetacellsByGroupsArgs$ident.group %||% 'seurat_clusters'}`) ",
        "and normalized with `NormalizeMetacells`."
    )
}

downstream_info <- c()
if (length(dmes) > 0) downstream_info <- c(downstream_info, glue("{length(dmes)} differential module expression analysis(es)"))
if (length(module_trait_corr) > 0) downstream_info <- c(downstream_info, glue("{length(module_trait_corr)} module-trait correlation analysis(es)"))
if (length(enrich) > 0) downstream_info <- c(downstream_info, glue("{length(enrich)} module enrichment analysis(es)"))
if (length(module_preservations) > 0) downstream_info <- c(downstream_info, glue("{length(module_preservations)} module preservation analysis(es)"))
if (length(tf_network) > 0) downstream_info <- c(downstream_info, "a TF regulatory network analysis")
if (length(motifs) > 0) downstream_info <- c(downstream_info, "a motif overlap analysis")

intro <- paste0(
    "hdWGCNA (high-dimensional Weighted Gene Co-expression Network Analysis) is a method for ",
    "constructing gene co-expression networks from single-cell RNA-seq data. The single-cell ",
    "expression data is aggregated into pseudobulk samples or metacells to alleviate the sparsity ",
    "and dropout issues, and the gene co-expression network is then constructed using the WGCNA ",
    "framework, where the genes are grouped into modules of co-expressed genes. The modules can ",
    "be further characterized by downstream analyses, such as differential module expression, ",
    "module-trait correlation, functional enrichment, module preservation, and TF regulatory ",
    "network analysis.\n\n",
    "In this analysis, the network was constructed as ", mode_str, " (WGCNA name: `", wgcna_name, "`). ",
    "The genes were selected with `SetupForWGCNA` (method: `", SetupForWGCNAArgs$gene_select %||% "fraction", "`",
    if (!is.null(SetupForWGCNAArgs$fraction)) paste0(", fraction: ", SetupForWGCNAArgs$fraction), "). ",
    metacell_info, " ",
    soft_power_info, " ",
    if (!is.null(mod_table)) {
        glue("{n_modules} non-grey modules were detected, containing {n_genes} genes in total.")
    } else {
        ""
    },
    " The TOM matrices were cached in the output directory (see the `cache` env) for reuse across reruns.",
    if (length(downstream_info) > 0) paste0(" The downstream analyses include: ", paste(downstream_info, collapse = ", "), ".") else "",
    "\n\nThis report contains the plots, tables, and descriptions of the entire analysis."
)

reporter$add2(
    list(name = "Overview", contents = list(
        list(kind = "descr", content = intro)
    )),
    hs = c("Introduction", "Overview"),
    ui = "tabs"
)
# the soft power and dendrogram plots (plots step) are merged into the
# Introduction section, so keep the Overview first within it
report <- reporter$report
report$Introduction <- c(
    report$Introduction["Overview"],
    report$Introduction[setdiff(names(report$Introduction), "Overview")]
)
reporter$report <- c(
    report["Introduction"],
    report[setdiff(names(report), "Introduction")]
)

log$info("Saving the Seurat object ...")
save_obj(srtobj, file.path(outdir, paste0(basename(outdir), ".qs")))

reporter$save(joboutdir)
