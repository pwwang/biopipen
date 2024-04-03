source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(rlang)
library(Seurat)
library(tidyseurat)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
joboutdir <- {{job.outdir | r}}  # nolint
mutaters <- {{envs.mutaters | r}}  # nolint
group.by <- {{envs["group-by"] | r}}  # nolint
ident.1 <- {{envs["ident-1"] | r}}  # nolint
ident.2 <- {{envs["ident-2"] | r}}  # nolint
each <- {{envs.each | r}}  # nolint
prefix_each <- {{envs.prefix_each | r}}  # nolint
subset <- {{envs.subset | r}}  # nolint
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

log_info("- Reading srtobj...")

srtobj <- readRDS(srtfile)
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    group.by = group.by,
    ident.1 = ident.1,
    ident.2 = ident.2,
    each = each,
    prefix_each = prefix_each,
    subset = subset,
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

expand_each <- function(name, case) {
    outcases <- list()
    if (is.null(case$each) || nchar(case$each) == 0) {
        if (is.null(case$section) || case$section == "DEFAULT") {
            outcases[[name]] <- case
        } else {
            outcases[[paste0(case$section, "::", name)]] <- case
        }
    } else {
        if (!is.null(case$section) && case$section != "DEFAULT") {
            log_warn("  Ignoring `section` in case `{name}` when `each` is set.")
            case$section <- NULL
        }
        if (is.null(case$subset)) {
            eachs <- srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            eachs <- srtobj@meta.data %>% dplyr::filter(!!!parse_exprs(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }
        for (each in eachs) {
            by <- make.names(paste0("..", name, "_", case$each,"_", each))
            srtobj@meta.data <<- srtobj@meta.data %>%
                mutate(!!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group.by),
                    NA
                ))

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

log_info("- Expanding cases...")
cases <- expand_cases(cases, defaults, expand_each)


ensure_sobj <- function(expr, allow_empty) {
    tryCatch({ expr }, error = function(e) {
        if (allow_empty) {
            log_warn("  Ignoring this case: {e$message}")
            return(NULL)
        } else {
            stop(e)
        }
    })
}


do_case <- function(name, case) {
    log_info("- Handling case: {name} ...")
    info <- casename_info(name, cases, outdir, create = TRUE)

    allow_empty = startsWith(case$group.by, "..")
    # prepare expression matrix
    log_info("  Preparing expression matrix...")
    sobj <- ensure_sobj({ srtobj %>% filter(!is.na(!!sym(case$group.by))) }, allow_empty)
    if (is.null(sobj)) { return() }

    if (!is.null(case$subset)) {
        sobj <- ensure_sobj({ sobj %>% filter(!!!parse_exprs(case$subset)) }, allow_empty)
        if (is.null(sobj)) { return() }
    }
    if (!is.null(case$ident.2)) {
        sobj <- ensure_sobj({ sobj %>% filter(!!sym(case$group.by) %in% c(case$ident.1, case$ident.2)) }, allow_empty)
        if (is.null(sobj)) { return() }
    }

    allclasses <- sobj@meta.data[, case$group.by, drop = TRUE]
    if (is.null(case$ident.2)) {
        case$ident.2 <- ".rest"
        allclasses[allclasses != case$ident.1] <- ".rest"
    }
    exprs <- GetAssayData(sobj, layer = "data")

    # get preranks
    log_info("  Getting preranks...")
    ranks <- prerank(exprs, case$ident.1, case$ident.2, allclasses, case$method)
    write.table(
        ranks,
        file.path(info$casedir, "fgsea.rank"),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )
    if (sum(is.na(ranks[, 2])) == nrow(ranks)) {
        if (length(allclasses) < 100) {
            log_warn("  Ignoring this case because all gene ranks are NA and there are <100 cells.")
            cat(
                paste0("Not enough cells (n = ", length(allclasses), ") to run fgsea."),
                file = file.path(info$casedir, "fgsea.log")
            )
            add_report(
                list(
                    kind = "error",
                    content = paste0("Not enough cells (n = ", length(allclasses), ") to run fgsea.")
                ),
                h1 = info$h1,
                h2 = info$h2
            )
            return()
        } else {
            stop(paste0(
                "All gene ranks are NA (# cells = ",
                length(allclasses),
                "). ",
                "It's probably due to high missing rate in the data. ",
                "You may want to try a different `envs$method` for pre-ranking."
            ))
        }
    }

    # run fgsea
    log_info("  Running fgsea...")
    case$rest$minSize <- case$minsize
    case$rest$maxSize <- case$maxsize
    case$rest$eps <- case$eps
    case$rest$nproc <- case$ncores
    runFGSEA(ranks, gmtfile, case$top, info$casedir, case$rest)

    add_report(
        list(kind = "fgsea", dir = info$casedir),
        h1 = info$h1,
        h2 = info$h2
    )
}

sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))

save_report(joboutdir)
