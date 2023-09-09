source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(rlang)
library(Seurat)
library(tidyseurat)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group.by <- {{envs["group-by"] | r}}  # nolint
ident.1 <- {{envs["ident-1"] | r}}  # nolint
ident.2 <- {{envs["ident-2"] | r}}  # nolint
each <- {{envs.each | r}}  # nolint
section <- {{envs.section | r}}  # nolint
gmtfile <- {{envs.gmtfile | r}}  # nolint
method <- {{envs.method | r}}  # nolint
top <- {{envs.top | r}}  # nolint
minsize <- {{envs.minsize | r}}  # nolint
maxsize <- {{envs.maxsize | r}}  # nolint
eps <- {{envs.eps | r}}  # nolint
ncores <- {{envs.ncores | r}}  # nolint
rest <- {{envs.rest | r: todot="-"}}  # nolint
cases <- {{envs.cases | r: todot="-"}}  # nolint

srtobj <- readRDS(srtfile)
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

expand_cases <- function() {
    # fill up cases with missing parameters
    if (is.null(cases) || length(cases) == 0) {
        filled_cases <- list(
            DEFAULT = list(
                group.by = group.by,
                ident.1 = ident.1,
                ident.2 = ident.2,
                each = each,
                section = section,
                gmtfile = gmtfile,
                method = method,
                top = top,
                minsize = minsize,
                maxsize = maxsize,
                eps = eps,
                ncores = ncores,
                rest = rest
            )
        )
    } else {
        filled_cases <- list()
        for (name in names(cases)) {
            case <- list_setdefault(
                cases[[name]],
                group.by = group.by,
                ident.1 = ident.1,
                ident.2 = ident.2,
                each = each,
                section = section,
                gmtfile = gmtfile,
                method = method,
                top = top,
                minsize = minsize,
                maxsize = maxsize,
                eps = eps,
                ncores = ncores,
                rest = rest
            )
            case$rest <- list_setdefault(case$rest, rest)
            filled_cases[[name]] <- case
        }
    }

    outcases <- list()
    # expand each
    for (name in names(filled_cases)) {
        case <- filled_cases[[name]]
        if (is.null(case$each) || nchar(case$each) == 0) {
            outcases[[paste0(case$section, ":", name)]] <- case
        } else {
            eachs <- srtobj@meta.data %>% pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
            for (each in eachs) {
                by <- make.names(paste0(".", name, "_", case$each,"_", each))
                srtobj@meta.data <<- srtobj@meta.data %>%
                    mutate(!!sym(by) := if_else(
                        !!sym(case$each) == each,
                        !!sym(case$group.by),
                        NA
                    ))

                key <- paste0(case$each, ":", each)
                if (name != "DEFAULT") {
                    key <- paste0(key, " - ", name)
                }
                outcases[[key]] <- case
                outcases[[key]]$group.by <- by
            }
        }
    }
    outcases
}

do_case <- function(name, case) {
    print(paste("- Processing case:", name, "..."))
    section_case <- unlist(strsplit(name, ":"))
    section <- section_case[1]
    casename <- paste(section_case[-1], collapse = ":")
    case_dir <- file.path(outdir, section, casename)
    dir.create(case_dir, showWarnings = FALSE, recursive = TRUE)

    # prepare expression matrix
    print("  Preparing expression matrix...")
    sobj <- srtobj %>% tidyseurat::filter(!is.na(!!sym(case$group.by)))
    if (!is.null(case$ident.2)) {
        sobj <- sobj %>% tidyseurat::filter(!!sym(case$group.by) %in% c(case$ident.1, case$ident.2))
    }

    allclasses <- sobj@meta.data[, case$group.by, drop = TRUE]
    if (is.null(case$ident.2)) {
        case$ident.2 <- ".rest"
        allclasses[allclasses != case$ident.1] <- ".rest"
    }
    exprs <- GetAssayData(sobj, slot = "data", assay = "RNA")

    # get preranks
    print("  Getting preranks...")
    ranks <- prerank(exprs, case$ident.1, case$ident.2, allclasses, case$method)
    write.table(
        ranks,
        file.path(case_dir, "fgsea.rank"),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )

    # run fgsea
    print("  Running fgsea...")
    case$rest$minSize <- case$minsize
    case$rest$maxSize <- case$maxsize
    case$rest$eps <- case$eps
    case$rest$nproc <- case$ncores
    runFGSEA(ranks, gmtfile, case$top, case_dir, case$rest)
}

cases <- expand_cases()
sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))
