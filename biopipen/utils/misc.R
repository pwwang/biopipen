# Misc utilities for R

.isBQuoted = function(x) {
    # Check if x is backtick-quoted
    nchar(x) >= 2 && x[1] == "`" && x[length(x)] == "`"
}

bQuote = function(x) {
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

list_setdefault = function(x, ...) {
    # Set the default value of a key in a list
    if (is.null(x)) {
        x <- list()
    }
    if (!is.list(x)) {
        stop("list_setdefault: list expected")
    }
    y = list(...)
    for (k in names(y)) {
        if (is.null(x[[k]])) {
            x[[k]] <- y[[k]]
        }
    }
    x
}

list_update = function(x, y) {
    # Update the value in x from y
    if (is.null(x)) {
        x <- list()
    }
    if (is.null(y)) {
        y <- list()
    }
    for (k in names(y)) {
        x[[k]] <- y[[k]]
    }
    x
}
