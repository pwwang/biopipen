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
