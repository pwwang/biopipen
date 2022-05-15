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
