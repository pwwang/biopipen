# Misc utilities for R
suppressPackageStartupMessages({
    library(logger)
    library(rlang)
    library(jsonlite)
})

.logger_layout <- layout_glue_generator(
    format = '{sprintf("%-7s", level)} [{format(time, "%Y-%m-%d %H:%M:%S")}] {msg}'
)
# print also debug messages, let pipen-poplog to filter
log_threshold(DEBUG)
log_layout(.logger_layout)
log_appender(appender_stdout)
tryCatch(log_errors(), error = function(e) {})

.isBQuoted <- function(x) {
    # Check if x is backtick-quoted
    nchar(x) >= 2 && startsWith(x, "`") && endsWith(x, "`")
}

bQuote <- function(x) {
    if (.isBQuoted(x)) {
        x
    } else {
        paste0("`", x, "`")
    }
}

.escape_regex <- function(x) {
    # Escape regex special characters
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

#' Slugify a string
#' Remember also update the one in gsea.R
#' @param x A string
#' @param non_alphanum_replace Replace non-alphanumeric characters
#' @param collapse_replace Collapse consecutive non-alphanumeric character replacements
#' @param tolower Convert to lowercase
#' @return A slugified string
slugify <- function(x, non_alphanum_replace="-", collapse_replace=TRUE, tolower=FALSE) {
    subs <- list(
        "š"="s", "œ"="oe", "ž"="z", "ß"="ss", "þ"="y", "à"="a", "á"="a", "â"="a",
        "ã"="a", "ä"="a", "å"="a", "æ"="ae", "ç"="c", "è"="e", "é"="e", "ê"="e",
        "ë"="e", "ì"="i", "í"="i", "î"="i", "ï"="i", "ð"="d", "ñ"="n", "ò"="o",
        "ó"="o", "ô"="o", "õ"="o", "ö"="o", "ø"="oe", "ù"="u", "ú"="u", "û"="u",
        "ü"="u", "ý"="y", "ÿ"="y", "ğ"="g", "ı"="i", "ĳ"="ij", "ľ"="l", "ň"="n",
        "ř"="r", "ş"="s", "ť"="t", "ų"="u", "ů"="u", "ý"="y", "ź"="z", "ż"="z",
        "ſ"="s", "α"="a", "β"="b", "γ"="g", "δ"="d", "ε"="e", "ζ"="z", "η"="h",
        "θ"="th", "ι"="i", "κ"="k", "λ"="l", "μ"="m", "ν"="n", "ξ"="x", "ο"="o",
        "π"="p", "ρ"="r", "σ"="s", "τ"="t", "υ"="u", "φ"="ph", "χ"="ch", "ψ"="ps",
        "ω"="o", "ά"="a", "έ"="e", "ή"="h", "ί"="i", "ό"="o", "ύ"="u", "ώ"="o",
        "ϐ"="b", "ϑ"="th", "ϒ"="y", "ϕ"="ph", "ϖ"="p", "Ϛ"="st", "ϛ"="st", "Ϝ"="f",
        "ϝ"="f", "Ϟ"="k", "ϟ"="k", "Ϡ"="k", "ϡ"="k", "ϰ"="k", "ϱ"="r", "ϲ"="s",
        "ϳ"="j", "ϴ"="th", "ϵ"="e", "϶"="p"
    )
    # replace latin and greek characters to the closest english character
    for (k in names(subs)) {
        x <- gsub(k, subs[[k]], x)
    }
    x <- gsub("[^[:alnum:]_]", non_alphanum_replace, x)
    if(collapse_replace) x <- gsub(
        paste0(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", non_alphanum_replace), "+"),
        non_alphanum_replace,
        x
    )
    if(tolower) x <- tolower(x)
    x
}

do_call <- function (what, args, quote = FALSE, envir = parent.frame())  {

  # source: Gmisc
  # author: Max Gordon <max@gforge.se>

  if (quote)
    args <- lapply(args, enquote)

  if (is.null(names(args)) ||
      is.data.frame(args)){
    argn <- args
    args <- list()
  }else{
    # Add all the named arguments
    argn <- lapply(names(args)[names(args) != ""], as.name)
    names(argn) <- names(args)[names(args) != ""]
    # Add the unnamed arguments
    argn <- c(argn, args[names(args) == ""])
    args <- args[names(args) != ""]
  }

  if (inherits(x = what, what = "character")){
    if(is.character(what)){
      fn <- strsplit(what, "[:]{2,3}")[[1]]
      what <- if(length(fn)==1) {
        get(fn[[1]], envir=envir, mode="function")
      } else {
        get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
      }
    }
    call <- as.call(c(list(what), argn))
  }else if (inherits(x = what, "function")){
    f_name <- deparse(substitute(what))
    call <- as.call(c(list(as.name(f_name)), argn))
    args[[f_name]] <- what
  }else if (inherits(x = what, what="name")){
    call <- as.call(c(list(what, argn)))
  }

  eval(call,
       envir = args,
       enclos = envir)

}

#' Save the plot into multiple formats
#'
#' @param plot The plot object
#' @param prefix The prefix of the file
#' @param formats The formats to save
#' @param bg The background color
#' @param devpars The device parameters
#' @export
save_plot <- function(plot, prefix, devpars = NULL, bg = "white", formats = c("png", "pdf")) {
    devpars <- devpars %||% list()
    devpars$res <- devpars$res %||% 100
    if (!is.null(attr(plot, "width"))) {
        devpars$width <- devpars$width %||% (attr(plot, "width") * devpars$res)
        devpars$height <- devpars$height %||% (attr(plot, "height") * devpars$res)
    } else {
        devpars$width <- devpars$width %||% 800
        devpars$height <- devpars$height %||% 600
    }

    old_dev <- grDevices::dev.cur()
    for (fmt in formats) {
        filename = paste0(prefix, ".", fmt)
        dev <- ggplot2:::plot_dev(fmt, filename, dpi = devpars$res)
        dim <- ggplot2:::plot_dim(c(devpars$width, devpars$height), units = "px", limitsize = FALSE, dpi = devpars$res)
        dev(filename = filename, width = dim[1], height = dim[2], bg = bg)
        print(plot)
        grDevices::dev.off()
    }
    on.exit(utils::capture.output({
        if (old_dev > 1) grDevices::dev.set(old_dev) # restore old device unless null device
    }))
}

#' Save the code to generate the data
#'
#' @param code The code
#' @param plot The plot object
#' @param setup The setup code to generate the plot
#' @param prefix The prefix of the file
#' @param ... Additional data frame to save
#'
#' @export
save_plotcode <- function(...) UseMethod("save_plotcode")

save_plotcode.character <- function(code, prefix, ..., envir = parent.frame()) {
    codedir <- paste0(prefix, ".code")
    dir.create(codedir, showWarnings = FALSE)
    codefile <- file.path(codedir, "plot.R")
    writeLines(code, codefile)
    save(..., file = file.path(codedir, "data.RData"), envir = envir)

    zip_file <- paste0(prefix, ".code.zip")
    zip::zip(zip_file, c("plot.R", "data.RData"), root = codedir)
    unlink(codedir, recursive = TRUE)
}

save_plotcode.ggplot <- function(plot, setup, prefix, ..., envir = parent.frame()) {
    if (is.null(plot$logs)) {
        stop("The plot object does not have logs, did you use gglogger?")
    }
    code <- plot$logs$gen_code(setup = setup)
    save_plotcode(code, prefix, ..., envir = envir)
}

#' Set the default value of a key in a list
#'
#' @param x A list
#' @param ... A list of key-value pairs
#' @return The updated list
#' @export
list_setdefault <- function(x, ...) {
    # Set the default value of a key in a list
    x <- x %||% list()

    stopifnot(is.list(x))
    y <- list(...)
    for (k in names(y)) {
        if (!k %in% names(x)) {
            # x[[k]] <- y[[k]]
            x <- c(x, y[k])
        }
    }
    x
}

#' Update a list with another list
#'
#' @param x A list
#' @param y A list
#' @param depth The depth to update, -1 means update all
#' @return The updated list
#' @export
list_update <- function(x, y, depth = -1L) {
    # Update the value in x from y
    x <- x %||% list()
    y <- y %||% list()

    for (k in names(y)) {
        if (is.null(y[[k]])) {
            x[[k]] <- NULL
            x <- c(x, y[k])
        } else if (is.list(x[[k]]) && is.list(y[[k]]) && depth != 0L) {
            x[[k]] <- list_update(x[[k]], y[[k]], depth - 1L)
        } else {
            x[[k]] <- y[[k]]
        }
    }
    x
}

#’ Biopipen palette
#'
#’ @param alpha The alpha value
#’ @return A palette function
#' @export
pal_biopipen <- function(alpha = 1) {
    if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
    colors <- c(
        "#ec3f3f", "#009e73", "#008ad8", "#cc79a7",
        "#e69f00", "#50cada", "#f0e442", "#a76ce7",
        "#ff864d", "#45e645", "#3699b5", "#ffdcda",
        "#d55e00", "#778ba6", "#c37b35", "#bc28ff"
    )
    colors <- scales::alpha(colors, alpha)
    function(n) {
        if (n <= length(colors)) {
            colors[1:n]
        } else {
            out_colors <- colors
            out_alpha <- 1.0
            while(length(out_colors) < n) {
                out_alpha <- out_alpha - 0.3
                out_colors <- c(out_colors, scales::alpha(colors, out_alpha))
            }
            out_colors[1:n]
        }
    }
}

scale_color_biopipen <- function(alpha = 1, ...) {
    ggplot2::discrete_scale("colour", "biopipen", pal_biopipen(alpha), ...)
}

scale_colour_biopipen <- scale_color_biopipen

scale_fill_biopipen <- function(alpha = 1, ...) {
    ggplot2::discrete_scale("fill", "biopipen", pal_biopipen(alpha), ...)
}

.report <- list(
    # h1 => list(
    #   h2 => list(
    #       h3#1 => list(ui1 => list(content11, content12)),
    #       h3#2 => list(ui2 => list(content21, content22))
    #   )
    # )
)

add_report <- function(..., h1, h2 = "#", h3 = "#", ui = "flat") {
    if (is.null(.report[[h1]])) {
        .report[[h1]] <<- list()
    }
    if (is.null(.report[[h1]][[h2]])) {
        .report[[h1]][[h2]] <<- list()
    }
    if (is.null(.report[[h1]][[h2]][[h3]])) {
        .report[[h1]][[h2]][[h3]] <<- list()
    }
    if (is.null(.report[[h1]][[h2]][[h3]][[ui]])) {
        .report[[h1]][[h2]][[h3]][[ui]] <<- list()
    }
    content = list(...)
    for (i in seq_along(content)) {
        .report[[h1]][[h2]][[h3]][[ui]] <<- c(
            .report[[h1]][[h2]][[h3]][[ui]],
            list(content[[i]])
        )
    }
}

save_report <- function(path, clear = TRUE) {
    if (dir.exists(path)) {
        path <- file.path(path, "report.json")
    }

    writeLines(toJSON(.report, pretty = TRUE, auto_unbox = TRUE), path)
    if (clear) {
        .report <<- list()
    }
}


# Escape html
html_escape <- function(text) {
    if (is.null(text)) { return("") }
    text <- gsub("&", "&amp;", text)
    text <- gsub("<", "&lt;", text)
    text <- gsub(">", "&gt;", text)
    text <- gsub("\"", "&quot;", text)
    text <- gsub("'", "&#039;", text)
    text
}

#' Expand the cases with default values
#' If a case has a key `each`, then it will be expanded by `expand_each`
#'
#' @param cases A list of cases
#' @param defaults A list of default values
#' @param expand_each A function to expand each case, if NULL, then the `each` key will be ignored.
#'   The function should take two arguments, `name` and `case`, and return a list of expanded cases.
#' @return A list of expanded cases
#' @export
expand_cases <- function(cases, defaults, expand_each = NULL) {
    if (is.null(cases) || length(cases) == 0) {
        filled_cases <- list(DEFAULT = defaults)
    } else {
        filled_cases <- list()
        for (name in names(cases)) {
            case <- list_update(defaults, cases[[name]], depth = 5L)
            filled_cases[[name]] <- case
        }
    }

    if (is.null(expand_each)) {
        return(filled_cases)
    }

    stopifnot(is.function(expand_each))

    outcases <- list()
    for (name in names(filled_cases)) {
        case <- filled_cases[[name]]
        each_cases <- expand_each(name, case)
        outcases <- c(outcases, each_cases)
    }
    outcases
}

#' Create information for a casename
#'
#' @param casename A casename
#' @param cases Used to check if there is only a single section in the cases
#' @param section_key The key to check, default is `section`
#' @param section The default section if no section if provided in the casename
#' @param sep The separator in casename to split section and casename
#' @param create Create the directory if not exists
#' @return A list of information, including `casedir`, `section`, `case`,
#'   `section_slug`, `case_slug` and the original `casename`.
#' @export
casename_info <- function(
    casename, cases, outdir,
    section_key = "section",
    section = NULL,
    sep = "::",
    case_type = c("dir", "prefix"),
    create = FALSE
) {
    section <- section %||% "DEFAULT"
    case_type <- match.arg(case_type)
    # CR_vs_PD_in_BL:seurat_clusters - IM IL1
    sec_case_names <- strsplit(casename, sep)[[1]]
    # seurat_clusters - IM IL1
    # In case we have more than one colon
    cname <- paste(sec_case_names[-1], collapse = "::")
    if (length(cname) == 0 || nchar(cname) == 0) {
        # no sep
        cname <- casename
    } else {
        section <- sec_case_names[1]
    }
    single_section <- length(unique(sapply(cases, function(x) x[[section_key]]))) == 1

    out <- list(
        casename = casename,
        section = section,
        case = cname,
        section_slug = slugify(section),
        case_slug = slugify(cname),
        h1 = ifelse(
            single_section && section == "DEFAULT",
            html_escape(cname),
            html_escape(ifelse(single_section, paste0(section, ": ", cname), section))
        ),
        h2 = ifelse(
            single_section && section == "DEFAULT",
            "#",
            ifelse(single_section, "#", html_escape(cname))
        )
    )

    if (case_type == "dir") {
        if (single_section && section == "DEFAULT") {
            out$casedir <- file.path(outdir, out$case_slug)
        } else {
            out$casedir <- file.path(outdir, out$section_slug, out$case_slug)
        }
        if (create) {
            dir.create(out$casedir, showWarnings = FALSE, recursive = TRUE)
        }
    } else {
        if (single_section && section == "DEFAULT") {
            out$caseprefix <- file.path(outdir, out$case_slug)
        } else {
            out$caseprefix <- file.path(outdir, out$section_slug, out$case_slug)
            if (create) {
                dir.create(file.path(outdir, out$section_slug), showWarnings = FALSE, recursive = TRUE)
            }
        }
    }
    out
}

run_command <- function(
    cmd,
    fg = FALSE,
    wait = TRUE,
    print_command = TRUE,
    print_command_handler = cat,
    ...
) {
    if (print_command) {
        print_command_handler("RUNNING COMMAND:\n")
        print_command_handler(paste0("  ", paste(cmd, collapse = " "), "\n\n"))
    }

    kwargs <- list(...)
    stdin <- kwargs$stdin %||% ""
    stdout <- kwargs$stdout %||% ""
    stderr <- kwargs$stderr %||% ""
    input <- kwargs$input %||% NULL
    k_env <- kwargs$env %||% list()
    env <- ""
    if (is.list(k_env)) {
        for (k in names(env)) { env <- paste0(env, k, "=", k_env[[k]], ";")}
    } else {
        env <- k_env
    }
    if (fg) {
        stdout <- ""
        stderr <- ""
    } else {
        if (stdout == "") { stdout <- FALSE }
    }

    command = cmd[1]
    args = cmd[-1]
    out <- system2(
        command,
        args = args,
        stdout = stdout,
        stderr = stderr,
        stdin = stdin,
        env = env,
        wait = wait,
        input = input
    )
    if (!isTRUE(stdout) && !isTRUE(stderr)) {
        if(out != 0) stop(sprintf("Command failed with exit code %s", out))
        if (!fg) { return(out) }
    } else {
        status <- attr(out, "status")
        if (is.integer(status) && status != 0) {
            stop(sprintf("Command failed: code (%s): %s", status, out))
        }
        return(out)
    }
}

#' Expand the dims usually used in single-cell analysis to specific dimensions
#'
#' @param dims The dimensions to expand
#' @return A vector of expanded dimensions
#' @export
#' @examples
#' expand_dims(NULL) # c(1, 2)
#' expand_dims(1:2) # c(1, 2)
#' expand_dims(1) # c(1)
#' expand_dims("1:2") # c(1, 2)
#' expand_dims("1") # c(1)
#' # dash works as the same as colon
#' expand_dims("1-3") # c(1, 2, 3)
#' expand_dims("1,3") # c(1, 3)
#' expand_dims("1,3:5") # c(1, 3, 4, 5)
#' expand_dims(c(1, "3:5", 7)) # c(1, 3, 4, 5, 7)
expand_dims <- function(dims, default = 1:2) {
    if (is.null(dims)) {
        return(default)
    }
    if (is.numeric(dims)) {
        return(dims)
    }
    dims <- unlist(strsplit(dims, ","))
    out <- c()
    for (d in dims) {
        if (grepl(":", d)) {
            d <- unlist(strsplit(d, ":"))
            d <- as.integer(d[1]):as.integer(d[2])
        } else if (grepl("-", d)) {
            d <- unlist(strsplit(d, "-"))
            d <- as.integer(d[1]):as.integer(d[2])
        } else {
            d <- as.integer(d)
        }
        out <- c(out, d)
    }
    out
}


#' Get plotthis function from plot_type
#'
#' @param plot_type The plot type
#' @param gglogger_register Register the plotthis function to gglogger
#' @param return_name Return the name of the function instead of the function
#' @return The plotthis function
#' @export
get_plotthis_fn <- function(plot_type, gglogger_register = TRUE, return_name = FALSE) {
    fn_name <- switch(plot_type,
        hist = "Histogram",
        histo = "Histogram",
        histogram = "Histogram",
        featuredim = "FeatureDimPlot",
        splitbar = "SplitBarPlot",
        enrichmap = "EnrichMap",
        enrichnet = "EnrichNetwork",
        enrichnetwork = "EnrichNetwork",
        gsea = "GSEAPlot",
        gseasummary = "GSEASummaryPlot",
        gseasum = "GSEASummaryPlot",
        heatmap = "Heatmap",
        network = "Network",
        pie = "PieChart",
        wordcloud = "WordCloudPlot",
        venn = "VennDiagram",
        paste0(tools::toTitleCase(plot_type), "Plot")
    )
    if (return_name) {
        return(fn_name)
    }
    fn <- tryCatch({
        utils::getFromNamespace(fn_name, "plotthis")
    }, error = function(e) {
        stop("Unknown plot type: ", plot_type)
    })

    if (gglogger_register) {
        gglogger::register(fn, fn_name)
    } else {
        fn
    }
}


#' Extract variables from a named list
#'
#' @param x A named list
#' @param ... The names of the variables
#' @param keep Keep the extracted variables in the list
#' @param env The environment to assign the extracted variables
#' @return The list with/ithout the extracted variables
#'
#' @export
extract_vars <- function(x, ..., keep = FALSE, env = parent.frame()) {
    stopifnot("extract_vars: 'x' must be a named list" = is.list(x) && !is.null(names(x)))
    vars <- list(...)
    if (is.null(names(vars))) {
        names(vars) <- unlist(vars)
    }
    for (i in seq_along(vars)) {
        if (nchar(names(vars)[i]) == 0) {
            names(vars)[i] <- vars[[i]]
        }
    }
    # list2env?
    for (n in names(vars)) {
        if (!n %in% names(x)) {
            stop(sprintf("Variable '%s' not found in the list", n))
        }
        assign(vars[[n]], x[[n]], envir = env)
        if (!isTRUE(keep)) {
            x[[n]] <- NULL
        }
    }

    x
}
