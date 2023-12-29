suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(tidyselect))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

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
#' @param each A column name (without quotes) in metadata to split the cells.
#'  Each comparison will be done for each value in this column.
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
#' @param debug Return the transformed data frame with counts, predicates, sum, and diff.
#' @param order The order of the returned ids. It could be `sum` or `diff`,
#'  which is the sum or diff of the `compare` between idents. Two kinds of
#'  modifiers can be added, including `desc` and `abs`. For example,
#'  `sum,desc` means the sum of `compare` between idents in descending order.
#'  Default is `diff,abs,desc`.
#'  It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
#'  ids will be in the same order as in `df`.
#' @param include_emerged Whether to include emerged clones for the expanded clones.
#'  Default is `FALSE`. It only works for `"expanded"`.
#' @param include_vanished Whether to include vanished clones for the collapsed clones.
#'  Default is `FALSE`. It only works for `"collapsed"`.
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
    each,
    uniq,
    order,
    debug
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
            "`compare` must be either a column name in df, or 'count'/'.n'. ",
            'Got "',
            compare_label,
            '"'
        ))
    }

    predicate <- function(ident_1, ident_2) {
        if (fun == "expanded") {
            ident_1 > ident_2 && ident_2 > 0
        } else if (fun == "expanded+") {
            ident_1 > ident_2
        } else if (fun == "collapsed") {
            ident_1 < ident_2 && ident_1 > 0
        } else if (fun == "collapsed+") {
            ident_1 < ident_2
        } else if (fun == "emerged") {
            ident_1 > 0 && ident_2 == 0
        } else if (fun == "vanished") {
            ident_1 == 0 && ident_2 > 0
        }
    }

    # subset the data frame
    trans <- df %>%
        dplyr::filter(!!subset) %>%
        drop_na(!!id) %>%
        # # remove NA values in group.by column
        # dplyr::filter(!is.na(!!group.by)) %>%
        # mark the group.by column (as .group) as ident_1 or ident_2 or NA
        mutate(
            .group = if_else(
                !!group.by == ident_1,
                "ident_1",
                if_else(ident_2 != "<NULL>" & !!group.by != ident_2, NA, "ident_2")
            )
        ) %>%
        # remove NA values in ..group column
        drop_na(.group)

    if (is_empty(attr(each, ".Environment"))) {
        if (as_label(each) == "NULL") {
            each <- NULL
        } else {
            each <- sym(as_name(each))
        }
    }
    if (is.null(each)) {
        trans <- trans %>% group_by(!!id, .group)
    } else {
        trans <- trans %>% group_by(!!each, !!id, .group)
    }

    if (compare_is_count) {
        trans <- trans %>% summarise(.n = n(), .groups = "drop")
    } else {
        trans <- trans %>% summarise(.n = first(!!compare), .groups = "drop")
    }

    trans <- trans %>% pivot_wider(names_from = .group, values_from = .n) %>%
        replace_na(list(ident_1 = 0, ident_2 = 0)) %>%
        rowwise() %>%
        # add the predicates, sums and diffs
        mutate(
            .predicate = predicate(ident_1, ident_2),
            .sum = ident_1 + ident_2,
            .diff = ident_1 - ident_2
        ) %>%
        ungroup() %>%
        arrange(!!order)

    if (debug) {
        return(trans)
    }

    uniq_ids <- trans %>% filter(.predicate) %>% pull(!!id) %>% as.vector() %>% unique()
    if (uniq) {
        return(uniq_ids)
    }

    df %>%
        mutate(
            .group = if_else(
                !!group.by == ident_1,
                "ident_1",
                if_else(ident_2 != "<NULL>" & !!group.by != ident_2, NA, "ident_2")
            ),
            .out = if_else(!!id %in% uniq_ids & !!subset & !is.na(.group), !!id, NA)
        ) %>%
        pull(.out)
}

#' @export
expanded <- function(
    df = .,
    group.by, # nolint
    idents,
    subset = TRUE,
    each = NULL,
    id = CDR3.aa,
    compare = .n,
    uniq = TRUE,
    debug = FALSE,
    order = desc(.sum),
    include_emerged = FALSE
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    fun = if (include_emerged) "expanded+" else "expanded"
    .size_compare(
        df = df,
        group.by = enquo(group.by),
        idents = idents,
        subset = enexpr(subset),
        id = enquo(id),
        compare = enquo(compare),
        fun = fun,
        each = tryCatch(enquo(each), error = function(e) NULL),
        uniq = uniq,
        order = enexpr(order),
        debug = debug
    )
}

#' @export
collapsed <- function(
    df = .,
    group.by, # nolint
    idents,
    subset = TRUE,
    each = NULL,
    id = CDR3.aa,
    compare = .n,
    uniq = TRUE,
    debug = FALSE,
    order = desc(.sum),
    include_vanished = FALSE
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    fun = if (include_vanished) "collapsed+" else "collapsed"
    .size_compare(
        df = df,
        group.by = enquo(group.by),
        idents = idents,
        subset = enexpr(subset),
        id = enquo(id),
        compare = enquo(compare),
        fun = fun,
        each = tryCatch(enquo(each), error = function(e) NULL),
        uniq = uniq,
        order = enexpr(order),
        debug = debug
    )
}

#' @export
emerged <- function(
    df = .,
    group.by, # nolint
    idents,
    subset = TRUE,
    each = NULL,
    id = CDR3.aa,
    compare = .n,
    uniq = TRUE,
    debug = FALSE,
    order = desc(.sum)
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df = df,
        group.by = enquo(group.by),
        idents = idents,
        subset = enexpr(subset),
        id = enquo(id),
        compare = enquo(compare),
        fun = "emerged",
        each = tryCatch(enquo(each), error = function(e) NULL),
        uniq = uniq,
        order = enexpr(order),
        debug = debug
    )
}

#' @export
vanished <- function(
    df = .,
    group.by, # nolint
    idents,
    subset = TRUE,
    each = NULL,
    id = CDR3.aa,
    compare = .n,
    uniq = TRUE,
    debug = FALSE,
    order = desc(.sum)
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }
    .size_compare(
        df = df,
        group.by = enquo(group.by),
        idents = idents,
        subset = enexpr(subset),
        id = enquo(id),
        compare = enquo(compare),
        fun = "vanished",
        each = tryCatch(enquo(each), error = function(e) NULL),
        uniq = uniq,
        order = enexpr(order),
        debug = debug
    )
}

#' Get paired entities from a data frame based on the other column
#'
#' @rdname Get paired entities
#' @param df The data frame. Use `.` if the function is called in a dplyr pipe.
#' @param id The column name in `df` for the ids to be returned in the
#'   final output
#' @param compare The column name in `df` to compare the values for each
#'   id in `id`.
#' @param idents The values in `compare` to compare. It could be either an
#'   an integer or a vector. If it is an integer, the number of values in
#'   `compare` must be the same as the integer for the `id` to be regarded
#'   as paired. If it is a vector, the values in `compare` must be the same
#'   as the values in `idents` for the `id` to be regarded as paired.
#' @param uniq Whether to return unique ids or not. Default is `TRUE`.
#'   If `FALSE`, you can mutate the meta data frame with the returned ids.
#'   Non-paired ids will be `NA`.
#' @return A vector of paired ids (in `id` column)
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
    df = .,
    id,
    compare,
    idents = 2,
    uniq = TRUE
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }

    id <- enquo(id)
    compare <- enquo(compare)
    if (is_empty(attr(id, ".Environment"))) {
        id <- sym(as_name(id))
    }
    if (is_empty(attr(compare, ".Environment"))) {
        compare <- sym(as_name(compare))
    }
    if (!as_name(id) %in% colnames(df)) {
        stop(paste0(
            '`id` must be a column name in df. Got "',
            as_name(id),
            '"'
        ))
    }
    if (!as_name(compare) %in% colnames(df)) {
        stop(paste0(
            '`compare` must be a column name in df. Got "',
            as_name(compare),
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
            add_count(!!id, name = "..count") %>%
            mutate(..paired = if_else(..count == idents, !!id, NA))
    } else {
        if (length(idents) <= 1) {
            stop(paste0(
                '`idents` must be a vector with length greater than 1. Got ',
                length(idents)
            ))
        }
        out <- df %>%
            group_by(!!id) %>%
            mutate(
                ..paired = if_else(
                    rep(setequal(!!compare, idents), n()),
                    !!id,
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

#' @export
#' @rdname Get top entities by size of group
#' @param df The data frame. Use `.` if the function is called in a dplyr pipe.
#' @param id The column name in `df` for the groups.
#' @param compare The column name in `df` to compare the values for each group.
#'   It could be either a numeric column or `.n` to compare the number of
#'   entities in each group. If a column is passed, the values in the column
#'   must be numeric and the same in each group. This won't be checked.
#' @param n The number of top entities to return. if `n` < 1, it will be
#'  regarded as the percentage of the total number of entities in each group
#'  (after subsetting or each applied).
#'  Specify 0 to return all entities.
#' @param subset An expression to subset the entities, will be passed to
#'   `dplyr::filter()`. Default is `TRUE` (no filtering).
#' @param with_ties Whether to return all entities with the same size as the
#'  last entity in the top list. Default is `FALSE`.
#' @param each A column name (without quotes) in metadata to split the cells.
#' @param debug Return the transformed data frame with counts and predicates
#' @param uniq Whether to return unique ids or not. Default is `TRUE`.
#'   If `FALSE`, you can mutate the meta data frame with the returned ids.
top <- function(
    df = .,
    id = CDR3.aa,
    n = 10,
    compare = .n,
    subset = TRUE,
    with_ties = FALSE,
    each = NULL,
    debug = FALSE,
    uniq = TRUE
) {
    lbl <- as_label(enquo(df))
    if (length(lbl) == 1 && lbl == ".") {
        df <- across(everything())
    }

    id <- enquo(id)
    compare <- enquo(compare)
    if (is.character(subset)) {
        subset <- parse_expr(subset)
    } else {
        subset <- enexpr(subset)
    }

    each <- tryCatch(enquo(each), error = function(e) NULL)
    if (is_empty(attr(id, ".Environment"))) {
        id <- sym(as_name(id))
    }
    if (is_empty(attr(compare, ".Environment"))) {
        compare <- sym(as_name(compare))
    }
    if (!as_name(id) %in% colnames(df)) {
        stop(paste0(
            '`id` must be a column name in df. Got "',
            as_name(id),
            '"'
        ))
    }
    if (!as_name(compare) %in% colnames(df) && as_name(compare) != '.n') {
        stop(paste0(
            '`compare` must be a column name in df. Got "',
            as_name(compare),
            '"'
        ))
    }
    if (is_empty(attr(each, ".Environment"))) {
        if (as_label(each) == "NULL") {
            each <- NULL
        } else {
            each <- sym(as_name(each))
        }
    }
    if (!is.null(each) && !as_name(each) %in% colnames(df)) {
        stop(paste0(
            '`each` must be a column name in df. Got "',
            as_name(each),
            '"'
        ))
    }

    subdf <- df %>% dplyr::filter(!!subset) %>% tidyr::drop_na(!!id)

    handle_one_each <- function(d) {
        if (!is.null(each)) {
            d <- d %>% group_by(!!each, !!id)
        } else {
            d <- d %>% group_by(!!id)
        }
        d <- d %>%
            dplyr::summarise(.n = dplyr::n(), .groups = "drop") %>%
            dplyr::arrange(dplyr::desc(!!compare))

        if (n > 0 && n < 1) {
            o <- d %>% dplyr::slice_max(prop = n, order_by = !!compare, with_ties = with_ties)
        } else if (n >= 1) {
            o <- d %>% dplyr::slice_max(n = n, order_by = !!compare, with_ties = with_ties)
        } else {
            o <- d
        }
        d %>% dplyr::mutate(.predicate = !!id %in% dplyr::pull(o, !!id))
    }

    if (is.null(each)) {
        out <- handle_one_each(subdf)
    } else {
        out <- subdf %>% dplyr::group_by(!!each) %>%
            dplyr::group_split() %>%
            purrr::map(handle_one_each) %>%
            dplyr::bind_rows()
    }

    if (isTRUE(debug)) {
        return(out)
    }

    uniq_ids <- out %>% dplyr::filter(.predicate) %>%
        dplyr::pull(!!id) %>% as.vector() %>% unique()
    if (isTRUE(uniq)) {
        return(uniq_ids)
    }

    df <- df %>% left_join(
        out,
        by = if(is.null(each)) as_name(id) else c(as_name(each), as_name(id)))

    df %>% dplyr::mutate(
        .out = if_else(.predicate & !!subset, !!id, NA)
    ) %>% dplyr::pull(.out)
}