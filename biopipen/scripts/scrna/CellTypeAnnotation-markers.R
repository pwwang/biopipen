# CellTypeAnnotation-markers.R — pure R functions, no Jinja2 template variables
# Source'd by CellTypeAnnotation.R
#
# Universal marker format: a long table (TSV/CSV/RDS/qs/qs2 data.frame) with
# columns `cell_type` and `gene` (required), plus optional `direction`
# (positive/negative), `weight`, `species`, `cancer`, `tissue`, and `level`.
# Column aliases are auto-detected. `file://path#col1,col2,...` selects
# columns by name (aliases allowed) or 1-based index. Files without
# `cell_type`/`gene` columns are returned as-is (native-format passthrough).

.MARKER_ALIASES <- list(
    cell_type = c("cell_type", "celltype", "type"),
    gene = c("gene", "marker", "gene_symbol"),
    direction = c("direction", "sign"),
    weight = c("weight"),
    species = c("species"),
    cancer = c("cancer"),
    tissue = c("tissue", "tissueType"),
    level = c("level")
)

canonicalize_marker_cols <- function(df) {
    if (!is.data.frame(df)) { return(df) }
    cn <- colnames(df)
    low_cn <- tolower(cn)
    for (canonical in names(.MARKER_ALIASES)) {
        if (canonical %in% cn) { next }
        hit <- which(low_cn %in% tolower(.MARKER_ALIASES[[canonical]]))[1]
        if (!is.na(hit)) { cn[hit] <- canonical }
    }
    colnames(df) <- cn
    df
}

is_marker_canonical <- function(df) {
    is.data.frame(df) &&
        "cell_type" %in% colnames(df) &&
        "gene" %in% colnames(df)
}

normalize_marker_direction <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[x %in% c("pos", "positive", "+")] <- "positive"
    x[x %in% c("neg", "negative", "-")] <- "negative"
    bad <- !x %in% c("positive", "negative")
    if (any(bad)) {
        stop(paste0(
            "Invalid `direction` value(s) in marker table: ",
            paste(unique(x[bad]), collapse = ", "),
            ". Accepted values: positive/negative (aliases: pos/neg/+/-)."
        ))
    }
    x
}

# Canonicalize columns and select a subset by name (exact, case-insensitive,
# alias) or 1-based index
apply_marker_cols <- function(df, cols, path) {
    if (!is.data.frame(df)) {
        stop(paste0("Cannot select columns from a non-table marker file: ", path))
    }
    df <- canonicalize_marker_cols(df)
    if (is.null(cols)) { return(df) }
    cn <- colnames(df)
    lookup <- setNames(cn, tolower(cn))
    for (canonical in names(.MARKER_ALIASES)) {
        if (canonical %in% cn) {
            for (alias in .MARKER_ALIASES[[canonical]]) {
                key <- tolower(alias)
                if (!key %in% names(lookup)) { lookup[[key]] <- canonical }
            }
        }
    }
    sel <- vapply(cols, function(col) {
        key <- tolower(col)
        if (key %in% names(lookup)) {
            lookup[[key]]
        } else if (grepl("^[0-9]+$", col)) {
            idx <- as.integer(col)
            if (idx < 1 || idx > length(cn)) {
                stop(paste0(
                    "Column index out of range in marker file: ", path, "#", col
                ))
            }
            cn[idx]
        } else {
            stop(paste0("Column not found in marker file: ", path, "#", col))
        }
    }, character(1))
    df[, sel, drop = FALSE]
}

# Read a marker file and return either a canonical marker data.frame
# (with `cell_type` and `gene` columns) or the raw object (a named list /
# matrix for RDS files, or the path string for xlsx files).
load_marker_table <- function(path, cols = NULL) {
    if (is.null(path)) { stop("The marker file path is not set") }

    if (startsWith(path, "file://")) {
        path <- sub("^file://", "", path)
    }
    if (is.null(cols) && grepl("#", path, fixed = TRUE)) {
        parts <- strsplit(path, "#", fixed = TRUE)[[1]]
        path <- parts[1]
        cols <- trimws(strsplit(parts[2], ",", fixed = TRUE)[[1]])
        cols <- cols[cols != ""]
        if (length(cols) == 0) { cols <- NULL }
    }

    if (!file.exists(path)) {
        stop(paste0("Marker file does not exist: ", path))
    }

    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) {
        return(path)  # native ScType xlsx passthrough
    }
    if (ext %in% c("rds", "qs", "qs2")) {
        obj <- read_obj(path)
        if (is.data.frame(obj)) {
            return(apply_marker_cols(obj, cols, path))
        }
        if (!is.null(cols)) {
            stop(paste0("Cannot select columns from a non-table marker file: ", path))
        }
        return(obj)
    }
    log$info("Loading marker table from {basename(path)} ...")
    df <- read.table(
        path,
        header = TRUE,
        sep = if (ext == "csv") "," else "\t",
        stringsAsFactors = FALSE
    )
    apply_marker_cols(df, cols, path)
}

# Convert a canonical marker table to the ScType format (one row per cell type)
markers_to_sctype_df <- function(df) {
    direction <- if ("direction" %in% colnames(df)) {
        normalize_marker_direction(df$direction)
    } else {
        rep("positive", nrow(df))
    }
    do.call(rbind, lapply(unique(df$cell_type), function(ct) {
        rows <- df$cell_type == ct
        pos <- sort(unique(df$gene[rows & direction == "positive"]))
        neg <- sort(unique(df$gene[rows & direction == "negative"]))
        tissue <- if ("tissue" %in% colnames(df)) {
            tv <- unique(df$tissue[rows])
            tv[!is.na(tv)][1] %||% NA
        } else {
            NA
        }
        level <- if ("level" %in% colnames(df)) {
            lv <- unique(df$level[rows])
            lv[!is.na(lv)][1] %||% 1
        } else {
            1
        }
        data.frame(
            tissueType = tissue,
            cellName = ct,
            geneSymbolmore1 = paste0(pos, collapse = ","),
            geneSymbolmore2 = paste0(neg, collapse = ","),
            Level = level,
            level = level,
            stringsAsFactors = FALSE
        )
    }))
}

# Convert a canonical marker table to the scSorter format
# Negative-direction markers become negative weights
markers_to_scsorter_df <- function(df) {
    anno <- data.frame(
        Type = df$cell_type,
        Marker = df$gene,
        stringsAsFactors = FALSE
    )
    has_dir <- "direction" %in% colnames(df)
    has_w <- "weight" %in% colnames(df)
    if (has_dir || has_w) {
        w <- if (has_w) as.numeric(df$weight) else rep(1, nrow(df))
        if (has_dir) {
            neg <- normalize_marker_direction(df$direction) == "negative"
            w[neg] <- -abs(w[neg])
            w[!neg] <- abs(w[!neg])
        }
        anno$Weight <- w
    }
    anno
}

# Keep only positive markers of a canonical marker table
filter_positive_markers <- function(df) {
    if ("direction" %in% colnames(df)) {
        df <- df[
            normalize_marker_direction(df$direction) == "positive",
            , drop = FALSE
        ]
    }
    df
}

# Convert a canonical marker table to a named list (cell type → genes)
# Only positive markers are kept (negative ones cannot be represented)
markers_to_named_list <- function(df) {
    df <- filter_positive_markers(df)
    split(as.character(df$gene), df$cell_type)
}

# Convert a canonical marker table to the scCATCH format
# Only positive markers are kept (negative ones cannot be represented)
markers_to_sccatch_df <- function(df) {
    df <- filter_positive_markers(df)
    colnames(df)[colnames(df) == "cell_type"] <- "celltype"
    df
}
