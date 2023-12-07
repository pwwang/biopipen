# Misc utilities for R
library(logger)
library(jsonlite)

.logger_layout <- layout_glue_generator(
    format = '{sprintf("%-7s", level)} [{format(time, "%Y-%m-%d %H:%M:%S")}] {msg}'
)
log_layout(.logger_layout)
log_appender(appender_stdout)
tryCatch(log_errors(), error = function(e) {})

.isBQuoted <- function(x) {
    # Check if x is backtick-quoted
    nchar(x) >= 2 && x[1] == "`" && x[length(x)] == "`"
}

bQuote <- function(x) {
    if (.isBQuoted(x)) {
        x
    } else {
        paste0("`", x, "`")
    }
}

slugify <- function(x, non_alphanum_replace="", space_replace="_", tolower=TRUE) {
  x <- gsub("[^[:alnum:] ]", non_alphanum_replace, x)
  x <- trimws(x)
  x <- gsub("[[:space:]]", space_replace, x)

  if(tolower) { x <- tolower(x) }

  return(x)
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

list_setdefault <- function(x, ...) {
    # Set the default value of a key in a list
    if (is.null(x)) {
        x <- list()
    }
    if (!is.list(x)) {
        stop("list_setdefault: list expected")
    }
    y <- list(...)
    for (k in names(y)) {
        if (!k %in% names(x)) {
            # x[[k]] <- y[[k]]
            x <- c(x, y[k])
        }
    }
    x
}

list_update <- function(x, y) {
    # Update the value in x from y
    if (is.null(x)) {
        x <- list()
    }
    if (is.null(y)) {
        y <- list()
    }
    for (k in names(y)) {
        if (is.null(y[[k]])) {
            x[[k]] <- NULL
            x <- c(x, y[k])
        } else {
            x[[k]] <- y[[k]]
        }
    }
    x
}

#’ Biopipen palette
#’ @param alpha Alpha value
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
    if (is.null(text)) {
        return("")
    }
    text = gsub("&", "&amp;", text)
    text = gsub("<", "&lt;", text)
    text = gsub(">", "&gt;", text)
    text = gsub("\"", "&quot;", text)
    text = gsub("'", "&#039;", text)
    text
}
