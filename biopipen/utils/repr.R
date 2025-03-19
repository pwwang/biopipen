library(rlang)

#' The string representation of an object
#' @param x An object
#' @param newline Whether to add newlines to the output for each element
#' @return The string representation
#' @export
repr <- function(x, newline = FALSE, ...) UseMethod("repr")

repr.default <- function(x, newline = FALSE, ...) {
    klass <- paste0(class(x), collapse = "/")
    fallback <- paste0("<", klass, ": ", deparse(substitute(x)), ">")

    tryCatch(
        x$.repr(newline, ...),
        error = function(e) {
            fallback
        }
    )
}

repr.numeric <- function(x, newline = FALSE, ...) {
    if (length(x) == 1) {
        as.character(x)
    } else if (!newline) {
        paste0("c(", paste(x, collapse = paste0(", ")), ")")
    } else {
        paste0(
            "c(\n",
            paste0(lapply(x, function(y) paste0("  ", y)), collapse = ",\n"),
            "\n)"
        )
    }
}

repr.character <- function(x, newline = FALSE, ...) {
    if (length(x) == 1) {
        paste0("\"", x, "\"")
    } else if (!newline) {
        paste0("c(", paste0(lapply(x, function(y) sQuote(y, q = FALSE)), collapse = ", "), ")")
    } else {
        paste0(
            "c(\n",
            paste0(lapply(x, function(y) paste0("  ", sQuote(y, q = FALSE))), collapse = ",\n"),
            "\n)"
        )
    }
}

repr.factor <- function(x, newline = FALSE, ...) {
    if (!newline) {
        paste0(
            "factor(", repr(as.character(x), newline, ...), ", levels = ", repr(levels(x), newline, ...), ")"
        )
    } else if (!newline) {
        paste0(
            "factor(\n",
            paste0("  ", repr(as.character(x), newline, ...), ",\n"),
            "  levels = ", repr(levels(x), newline, ...), "\n)"
        )
    } else {
        paste0(
            "factor(\n",
            paste0("  ", repr(as.character(x), newline, ...), ",\n"),
            "  levels = ", repr(levels(x), newline, ...), "\n)"
        )
    }
}

repr.logical <- function(x, newline = FALSE, ...) {
    if (length(x) == 1) {
        if (x) "TRUE" else "FALSE"
    } else if (!newline) {
        paste0("c(", paste0(x, collapse = ","), ")")
    } else {
        paste0(
            "c(\n",
            paste0(lapply(x, function(y) paste0("  ", y)), collapse = ",\n"),
            "\n)"
        )
    }
}

repr.list <- function(x, newline = FALSE, ...) {
  start <- if (newline) "list(\n" else "list("
  end <- if (newline) "\n)" else ")"
  sep <- if (newline) ",\n" else ", "
  prefix <- if (newline) "  " else ""
  if (length(names(x)) > 0) {
    paste0(
        start,
        paste0(
            lapply(seq_along(x), function(i) {
                name <- names(x)[i]
                if (identical(name, "")) {
                  paste0(prefix, repr(x[[i]]))
                } else {
                  paste0(prefix, bQuote(name), " = ", repr(x[[name]]))
                }
            }), collapse = sep
        ),
        end
    )
  } else {
    paste0(
        start, paste0(lapply(x, repr, newline = newline, ...), collapse = sep), end
    )
  }
}

repr.NULL <- function(x, newline = FALSE, ...) {
    "NULL"
}

repr.formula <- function(x, newline = FALSE, ...) {
  deparse(x)
}

repr.data.frame <- function(x, newline = FALSE, ...) {
  paste0(
      "data.frame(\n",
      paste0(
          lapply(names(x), function(name) {
              paste0("  ", bQuote(name), " = ", repr(x[[name]], newline = newline, ...))
          }), collapse = ",\n"
      ),
      "\n)"
  )
}

repr.environment <- function(x, newline = FALSE, ...) {
    if (!is_environment(x)) {
        # in case .GlobalEnv is dispatched here
        NextMethod()
    } else {
        nl <- if (newline) "\n" else " "
        prefix <- if (newline) "  " else ""
        paste0(
            "rlang::env(", nl, paste0(
                lapply(ls(x), function(name) {
                    paste0(prefix, bQuote(name), " = ", repr(get(name, x), newline = newline, ...))
                }), collapse = paste0(",", nl)
            ), nl, ")"
        )
    }
}
