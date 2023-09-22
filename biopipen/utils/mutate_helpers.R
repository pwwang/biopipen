suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tidyselect))
suppressPackageStartupMessages(library(dplyr))

#' Get expanded, collapsed, emerged or vanished clones from a meta data frame
#'
#' @rdname Get expanded, collapsed, emerged or vanished clones
#'
#' @param df The meta data frame
#' @param group.by The column name (without quotes) in metadata to group the
#'  cells.
#' @param idents The groups of cells to compare (values in `group-by` column).
#'  Either length 1 (`ident_1`) or length 2 (`ident_1` and `ident_2`).
#'  If length 1, the rest of the cells with non-NA values in `group.by` will
#'  be used as `ident_2`.
#' @param subset An expression to subset the cells, will be passed to
#'  `dplyr::filter()`. Default is `TRUE` (no filtering).
#' @param id The column name (without quotes) in metadata for the
#'  group ids (i.e. `CDR3.aa`)
#' @param compare Either a (numeric) column name (i.e. `Clones`, without quotes)
#'  in metadata to compare between groups, or `.n` to compare the
#'  number of cells in each group.
#' @param fun The way to compare between groups. Either `"expanded"`,
#'  `"collapsed"`, `"emerged"` or `"vanished"`.
#' @param uniq Whether to return unique ids or not. Default is `TRUE`.
#'  If `FALSE`, you can mutate the meta data frame with the returned ids.
#'  For example, `df %>% mutate(expanded = expanded(...))`.
#' @param order The order of the returned ids. It could be `sum` or `diff`,
#'  which is the sum or diff of the `compare` between idents. Two kinds of
#'  modifiers can be added, including `desc` and `abs`. For example,
#'  `sum,desc` means the sum of `compare` between idents in descending order.
#'  Default is `diff,abs,desc`.
#'  It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
#'  ids will be in the same order as in `df`.
#'
#' @return A vector of expanded or collapsed clones (in `id` column)
#'  If uniq is `FALSE`, the vector will be the same length as `df`.
#'
#' @examples
#' # Get expanded clones
#' df <- tibble(
#'  Clones = c(10, 8, 1, 5, 9, 2, 3, 7, 6, 4, 9, 9),
#'  Source = c(
#'      "Tumor", "Normal", "Normal", "Normal", "Tumor", "Tumor",
#'      "Tumor", "Normal", "Normal", "Normal", NA, "X"
#'  ),
#'  CDR3.aa = c("A", "C", "B", "E", "D", "E", "E", "B", "B", "B", "A", "A")
#' )
#'
#' expanded(df, Source, c("Tumor", "Normal"))
#' # The transformed data frame looks like this:
#   CDR3.aa ..predicate ..sum ..diff
#   <chr>   <lgl>       <dbl>  <dbl>
# 1 A       TRUE           10     10
# 2 B       FALSE           1     -1
# 3 C       FALSE           8     -8
# 4 D       TRUE            9      9
# 5 E       FALSE           7     -3
#'
#' # [1] "A" "D"
#'
#' # Get collapsed clones
#' collapsed(df, Source, c("Tumor", "Normal"))
#' # [1] "B" "C" "E"
#'
#' # Get emerged clones
#' emerged(df, Source, c("Tumor", "Normal"))
#' # [1] "A" "D"
#'
#' # Get vanished clones
#' vanished(df, Source, c("Tumor", "Normal"))
#' # [1] "B" "C"
.size_compare <- function(
    df,
    group.by, # nolint
    idents,
    subset,
    id,
    compare,
    fun,
    uniq,
    order
) {
    if (length(idents) == 1) {
        ident_1 <- idents[1]
        ident_2 <- NULL
    } else if (length(idents) == 2) {
        ident_1 <- idents[1]
        ident_2 <- idents[2]
    } else {
        stop("idents must be length 1 or 2")
    }
    if (is.null(ident_2)) ident_2 <- "<NULL>"

    if (is_empty(attr(group.by, ".Environment"))) {
        # Works if a (quoted) string passed
        group.by <- sym(as_name(group.by))
    }
    if (is_empty(attr(id, ".Environment"))) {
        id <- sym(as_name(id))
    }
    if (is_empty(attr(compare, ".Environment"))) {
        compare <- sym(as_name(compare))
    }
    compare_label <- as_name(compare)
    compare_is_count <- compare_label == '.n'

    if (!as_name(group.by) %in% colnames(df)) {
        stop(paste0(
            '`group.by` must be a column name in df. Got "',
            as_name(group.by),
            '"'
        ))
    }

    if (!compare_is_count && !compare_label %in% colnames(df)) {
        stop(paste0(
            "`compare` must be either a column name in df, or 'count'/'n'. ",
            'Got "',
            compare_label,
            '"'
        ))
    }

    predicate <- function(comp) {
        if (fun == "expanded") {
            comp[1] > comp[2]
        } else if (fun == "collapsed") {
            comp[1] < comp[2]
        } else if (fun == "emerged") {
            comp[1] > 0 && comp[2] == 0
        } else if (fun == "vanished") {
            comp[1] == 0 && comp[2] > 0
        }
    }

    # subset the data frame
    trans <- df %>% dplyr::filter(!!subset) %>%
        # remove NA values in group.by column
        dplyr::filter(!is.na(!!group.by)) %>%
        # mark the group.by column (as ..group) as ident_1 or ident_2 or NA
        mutate(
            ..group = if_else(
                !!group.by == ident_1,
                "ident_1",
                if_else(ident_2 != "<NULL>" & !!group.by != ident_2, NA, "ident_2")
            )
        ) %>%
        # remove NA values in ..group column
        dplyr::filter(!is.na(..group)) %>%
        # for each clone and group (ident_1 and ident_2)
        group_by(!!id, ..group) %>%
        # summarise the number of cells in each clone and group
        # so that we can compare between groups later
        summarise(
            ..compare = ifelse(compare_is_count, n(), first(!!compare)),
            .groups = "drop"
        ) %>%
        # for each clone, either compare Clones or ..count between groups
        # (ident_1 and ident_2)
        group_by(!!id) %>%
        # add missing group (either ident_1 or ident_2)
        group_modify(function(d, ...) {
            if (nrow(d) == 1) {
                d <- d %>% add_row(
                    ..group = ifelse(
                        d$..group == "ident_1", "ident_2", "ident_1"
                    ),
                    ..compare = 0
                )
            }
            d
        }) %>%
        # make sure ident_1 and ident_2 are in order
        arrange(..group, .by_group = TRUE) %>%
        # add the predicates, sums and diffs
        summarise(
            ..predicate = predicate(..compare),
            ..sum = sum(..compare),
            ..diff = ..compare[1] - ..compare[2]
        ) %>%
        # filter the clones
        dplyr::filter(..predicate)

    order_sum <- grepl("sum", order)
    order_diff <- grepl("diff", order)
    order_desc <- grepl("desc", order)
    order_abs <- grepl("abs", order)
    if (order_sum && !order_desc) {
        out <- trans %>% arrange(..sum) %>% pull(!!id)
    } else if (order_sum) {
        out <- trans %>% arrange(desc(..sum)) %>% pull(!!id)
    } else if (order_diff && !order_desc && !order_abs) {
        out <- trans %>% arrange(..diff) %>% pull(!!id)
    } else if (order_diff && !order_desc && order_abs) {
        out <- trans %>% arrange(abs(..diff)) %>% pull(!!id)
    } else if (order_diff && order_desc && !order_abs) {
        out <- trans %>% arrange(desc(..diff)) %>% pull(!!id)
    } else if (order_diff && order_desc && order_abs) {
        out <- trans %>% arrange(desc(abs(..diff))) %>% pull(!!id)
    } else {
        out <- trans %>% pull(!!id)
    }

    if (uniq) { return(out) }

    df %>% mutate(..out = if_else(!!id %in% out, !!id, NA)) %>% pull(..out)
}

#' @export
expanded <- function(
    df,
    group.by, # nolint
    idents,
    subset = TRUE,
    id = CDR3.aa,
    compare = Clones,
    uniq = TRUE,
    order = "diff+desc"
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df,
        enquo(group.by),
        idents,
        enquo(subset),
        enquo(id),
        enquo(compare),
        "expanded",
        uniq = uniq,
        order = order
    )
}

#' @export
collapsed <- function(
    df,
    group.by, # nolint
    idents,
    subset = TRUE,
    id = CDR3.aa,
    compare = Clones,
    uniq = TRUE,
    order = "diff+desc"
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df,
        enquo(group.by),
        idents,
        enquo(subset),
        enquo(id),
        enquo(compare),
        "collapsed",
        uniq = uniq,
        order = order
    )
}

#' @export
emerged <- function(
    df,
    group.by, # nolint
    idents,
    subset = TRUE,
    id = CDR3.aa,
    compare = Clones,
    uniq = TRUE,
    order = "diff+desc"
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df,
        enquo(group.by),
        idents,
        enquo(subset),
        enquo(id),
        enquo(compare),
        "emerged",
        uniq = uniq,
        order = order
    )
}

#' @export
vanished <- function(
    df,
    group.by, # nolint
    idents,
    subset = TRUE,
    id = CDR3.aa,
    compare = Clones,
    uniq = TRUE,
    order = "diff+desc"
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df,
        enquo(group.by),
        idents,
        enquo(subset),
        enquo(id),
        enquo(compare),
        "vanished",
        uniq = uniq,
        order = order
    )
}

#' Get paired entities from a data frame based on the other column
#'
#' @rdname Get paired entities
#' @param df The data frame. Use `.` if the function is called in a dplyr pipe.
#' @param id_col The column name in `df` for the ids to be returned in the
#'   final output
#' @param compare_col The column name in `df` to compare the values for each
#'   id in `id_col`.
#' @param idents The values in `compare_col` to compare. It could be either an
#'   an integer or a vector. If it is an integer, the number of values in
#'   `compare_col` must be the same as the integer for the `id` to be regarded
#'   as paired. If it is a vector, the values in `compare_col` must be the same
#'   as the values in `idents` for the `id` to be regarded as paired.
#' @param uniq Whether to return unique ids or not. Default is `TRUE`.
#'   If `FALSE`, you can mutate the meta data frame with the returned ids.
#'   Non-paired ids will be `NA`.
#' @return A vector of paired ids (in `id_col` column)
#' @examples
#' df <- tibble(
#'   id = c("A", "A", "B", "B", "C", "C", "D", "D"),
#'   compare = c(1, 2, 1, 1, 1, 2, 1, 2)
#' )
#' paired(df, id, compare, 2)
#' # [1] "A" "B" "C" "D"
#' paired(df, id, compare, c(1, 2))
#' # [1] "A" "C" "D"
#' paired(df, id, compare, c(1, 2), uniq = FALSE)
#' # [1] "A" "A" NA NA "C" "C" "D" "D"
#'
paired <- function(
    df,
    id_col,
    compare_col,
    idents = 2,
    uniq = TRUE
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }

    id_col <- enquo(id_col)
    compare_col <- enquo(compare_col)
    if (is_empty(attr(id_col, ".Environment"))) {
        id_col <- sym(as_name(id_col))
    }
    if (is_empty(attr(compare_col, ".Environment"))) {
        compare_col <- sym(as_name(compare_col))
    }
    if (!as_name(id_col) %in% colnames(df)) {
        stop(paste0(
            '`id_col` must be a column name in df. Got "',
            as_name(id_col),
            '"'
        ))
    }
    if (!as_name(compare_col) %in% colnames(df)) {
        stop(paste0(
            '`compare_col` must be a column name in df. Got "',
            as_name(compare_col),
            '"'
        ))
    }

    if (is.numeric(idents) && length(idents) == 1) {
        if (idents <= 1) {
            stop(paste0(
                '`idents` must be greater than 1. Got ',
                idents
            ))
        }
        out <- df %>%
            add_count(!!id_col, name = "..count") %>%
            mutate(..paired = if_else(..count == idents, !!id_col, NA))
    } else {
        if (length(idents) <= 1) {
            stop(paste0(
                '`idents` must be a vector with length greater than 1. Got ',
                length(idents)
            ))
        }
        out <- df %>%
            group_by(!!id_col) %>%
            mutate(
                ..paired = if_else(
                    rep(setequal(!!compare_col, idents), n()),
                    !!id_col,
                    NA
                )
            ) %>%
            ungroup()
    }

    out <- out %>% pull(..paired)
    if (uniq) {
        return(out %>% na.omit() %>% unique() %>% as.vector())
    } else {
        return(out)
    }
}
