source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(rlang)
library(Seurat)
library(tidyseurat)
library(slugify)

srtfile <- {{in.srtobj | r}}  # nolint
outdir <- {{out.outdir | r}}  # nolint
joboutdir <- {{job.outdir | r}}  # nolint
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

log_info("Reading srtobj...")

srtobj <- readRDS(srtfile)
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

single_section <- TRUE

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
    sections <- c()
    # expand each
    for (name in names(filled_cases)) {
        case <- filled_cases[[name]]
        if (is.null(case$each) || nchar(case$each) == 0) {
            sections <- c(sections, case$section)
            outcases[[paste0(case$section, ":", name)]] <- case
        } else {
            eachs <- srtobj@meta.data %>% pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
            sections <- c(sections, case$each)
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
    single_section <- length(unique(sections)) == 1
    outcases
}

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

do_case <- function(name, case) {
    log_info("- Doing case: {name} ...")
    info <- casename_info(name, create = TRUE)

    # prepare expression matrix
    log_info("  Preparing expression matrix...")
    sobj <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
    if (!is.null(case$ident.2)) {
        sobj <- sobj %>% filter(!!sym(case$group.by) %in% c(case$ident.1, case$ident.2))
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
                h1 = ifelse(
                    info$section == "DEFAULT",
                    info$case,
                    ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
                ),
                h2 = ifelse(
                    info$section == "DEFAULT",
                    "#",
                    ifelse(single_section, "#", info$case)
                )
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
        h1 = ifelse(
            info$section == "DEFAULT",
            info$case,
            ifelse(single_section, paste0(info$section, " - ", info$case), info$section)
        ),
        h2 = ifelse(
            info$section == "DEFAULT",
            "#",
            ifelse(single_section, "#", info$case)
        )
    )
}

cases <- expand_cases()
sapply(sort(names(cases)), function(name) do_case(name, cases[[name]]))

save_report(joboutdir)
