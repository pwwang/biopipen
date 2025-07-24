library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(hash)
library(glmnet)
library(broom.mixed)
library(stringr)
library(plotthis)
library(biopipen.utils)

scrfile <- {{in.scrfile | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
group_name <- {{envs.group | r}}
comparison <- {{envs.comparison | r}}
target <- {{envs.target | r}}
each_cols <- {{envs.each | r}}

log <- get_logger()
reporter <- get_reporter()

if (is.null(group_name) || is.null(comparison)) {
    stop("envs.group and envs.comparison must be specified")
}

if (length(comparison) != 2) {
    stop("envs.comparison must have exactly two elements or keys, representing the two groups to compare")
}

if (!is.list(comparison)) {
    comparison <- stats::setNames(as.list(comparison), comparison)
}

target <- target %||% names(comparison)[1]
if (!(target %in% names(comparison))) {
    stop(paste0("Target group '", target, "' not found in the comparison groups."))
}

if (is.character(each_cols) && length(each_cols) == 1) {
    each_cols = trimws(strsplit(each_cols, ",")[[1]])
}

### Helpers

get_positions <- function(seq_length, mid=FALSE){
  poss = c("104", "105", "106", "107", "108", "109", "110", "111", "111.1", "112.1", "112", "113", "114", "115", "116", "117", "118")
  if (mid){
    poss = poss[!(poss %in% c("104", "105", "106", "107", "113", "114", "115", "116", "117", "118"))]
  }
  if (seq_length==12){
    poss = poss[!(poss %in% c("110", "111", "111.1", "112.1", "112"))]
  }
  if (seq_length==13){
    poss = poss[!(poss %in% c("111", "111.1", "112.1","112"))]
  }
  if (seq_length==14){
    poss = poss[!(poss %in% c("111.1", "112.1", "111"))]
  }
  if (seq_length==15){
    poss = poss[!(poss %in% c("111.1", "112.1"))]
  }
  if (seq_length==16){
    poss = poss[!(poss %in% c("111.1"))]
  }
  return(poss)
}

get_feat_score <- function(x, amap){
  sum = 0
  for (i in 1:nchar(x)){
    sum = sum + amap[[substr(x, i, i)]]
  }
  return(sum/nchar(x))
}

add_percentAA <- function(data){
    data$midseq <- sapply(data$sequence, function(x) substr(x, 5, nchar(x)-6))
    for (i in 1:length(AA)){
        data$new <- sapply(data$midseq, function(x) 100*str_count(x, AA[i])/nchar(x))
        colnames(data)[ncol(data)] <- paste("perc_mid", AA[i], sep="_")
    }
    return(data)
}

add_positionalAA <- function(data, exc_VJ=FALSE){
  data$Vmotif = sapply(data$sequence, function(x) substr(x, 1, 3))
  poss = get_positions(17)
  for (i in 1:6){
    if (!(exc_VJ & i<4)){
      data$new <- sapply(data$sequence, function(x) substr(x, i, i))
      colnames(data)[ncol(data)] <- paste("p", poss[i], sep="")
    }
  }
  data$p110 <- sapply(data$sequence, function(x) ifelse(nchar(x)>12, substr(x, 7, 7), "*"))
  data$p111 <- sapply(data$sequence, function(x) ifelse(nchar(x)>14, substr(x, 8, 8), "*"))
  data$'p111.1' <- sapply(data$sequence, function(x) ifelse(nchar(x)>16, substr(x, 9, 9), "*"))
  data$'p112.1' <- sapply(data$sequence, function(x) ifelse(nchar(x)>15, substr(x, nchar(x)-7, nchar(x)-7), "*"))
  data$'p112' <- sapply(data$sequence, function(x) ifelse(nchar(x)>13, substr(x, nchar(x)-6, nchar(x)-6), "*"))
  for (i in 5:0){
    if (!(exc_VJ & i<5)){
      data$new <- sapply(data$sequence, function(x) substr(x, nchar(x)-i, nchar(x)-i))
      colnames(data)[ncol(data)] <- paste("p", poss[length(poss)-i], sep="")
    }
  }
  data$Jmotif = sapply(data$sequence, function(x) substr(x, nchar(x)-4, nchar(x)))
  jmotifs = c('AEAFF', 'AGYTF', 'DEQFF', 'DEQYF', 'DGYTF', 'DTQYF', 'EGYTF', 'EKLFF', 'ETQYF', 'GEAFF', 'GELFF', 'GEQFF', 'GEQYF', 'GGYTF', 'GKLFF', 'GTQYF', 'HEQFF', 'HEQYF', 'HGYTF', 'KTQYF', 'LGYTF', 'NEQFF', 'NEQYF', 'NGYTF', 'NIQYF', 'NTIYF', 'NTQYF', 'NVLTF', 'other', 'PEAFF', 'QPQHF', 'REQYF', 'RGYTF', 'SEAFF', 'SEQFF', 'SEQYF', 'SGYTF', 'SPLHF', 'VGYTF', 'YEQFF', 'YEQYF', 'YGYTF')
  data$Jmotif[!(data$Jmotif %in% jmotifs)] = "other"

  return(data)
}

create_hashmap <- function(key_vector, val_vector){
  h = hash()
  for (i in 1:length(key_vector)){
    h[[key_vector[i]]] = val_vector[i]
  }
  return(h)
}

# Configurable?
CDR3_MINLEN = 12
CDR3_MAXLEN = 17
AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
RF = tribble(
    ~AA,   ~pI,  ~IFH, ~volume,
    "A",  6.00, -0.17,   88.6,
    "R", 10.76, -0.81,  173.4,
    "N",  5.41, -0.42,  114.1,
    "D",  2.77, -1.23,  111.1,
    "C",  5.07,  0.24,  108.5,
    "Q",  5.65, -0.58,  143.8,
    "E",  3.22, -2.02,  138.4,
    "G",  5.97, -0.01,   60.1,
    "H",  7.59, -0.96,  153.2,
    "I",  6.02,  0.31,  166.7,
    "L",  5.98,  0.56,  166.7,
    "K",  9.74, -0.99,  168.6,
    "M",  5.74,  0.23,  162.9,
    "F",  5.48,  1.13,  189.9,
    "P",  6.30, -0.45,  112.7,
    "S",  5.68, -0.13,   89.0,
    "T",  5.60, -0.14,  116.1,
    "W",  5.89,  1.85,  227.8,
    "Y",  5.66,  0.94,  193.6,
    "V",  5.96, -0.07,  140.0
)
AA_FEATURES <- colnames(RF)[2:ncol(RF)]
AA_MAPS <- list()
for (i in 1:3){
  AA_MAPS[[i]] <- create_hashmap(as.character(RF$AA), as.vector(RF[,(i+1),drop=TRUE]))
}

log$info("Loading data from input file")
mdata <- read_obj(scrfile)@meta.data

if (!group_name %in% colnames(mdata)) {
    stop(paste0("Group name '", group_name, "' not found in the data."))
}

# check if valuess of comparison is in the group_name column
if (!all(unlist(comparison) %in% as.character(mdata[[group_name]]))) {
    stop(paste0("Some values in comparison are not found in the group_name column: ",
                paste(setdiff(unlist(comparison), mdata[[group_name]]), collapse = ", ")))
}

# add a new column with the keys of comparison, when their values are in the group_name column
mdata$.Group <- sapply(as.character(mdata[[group_name]]), function(x) {
    for (key in names(comparison)) {
        if (x %in% comparison[[key]]) {
            return(key)
        }
    }
    return(NA)
})
mdata <- mdata %>%
    separate(CTaa, into = c(NA, "sequence"), sep = "_", remove = FALSE) %>%
    separate(CTgene, into = c(NA, "vjgene"), sep = "_", remove = FALSE) %>%
    separate(vjgene, into = c("vgene", NA, "jgene", NA), sep = "\\.", remove = FALSE) %>%
    mutate(length = nchar(sequence))

# Statistics about the cell numbers with groups avaiable in metadata
# !!group_name, TotalCells, AvailCells, AvailCellsPct
log$info("Calculating statistics")
if (is.null(each_cols)) {
    stats = mdata %>%
        # group by group_name
        group_by(.Group) %>%
        summarise(
            TotalCells = nrow(mdata),
            CellsPerGroup = n(),
            AvailCellsPerGroup = sum(length >= CDR3_MINLEN & length <= CDR3_MAXLEN),
            # Percentage with % in character
            AvailCellsPct = paste0(round(AvailCellsPerGroup / CellsPerGroup * 100, 2), "%"),
            .groups = "drop"
        )
} else {
    stats = mdata %>%
        unite(".Subset", all_of(each_cols), sep = "_", remove = FALSE) %>%
        group_by(.Subset) %>%
        group_map(function(df, .y) {
            df %>%
                group_by(.Group) %>%
                summarise(
                    .Subset = .y$.Subset[1],
                    AllCells = nrow(mdata),
                    TotalCells = nrow(df),
                    CellsPerGroup = n(),
                    AvailCellsPerGroup = sum(length >= CDR3_MINLEN & length <= CDR3_MAXLEN),
                    # Percentage with % in character
                    AvailCellsPct = paste0(round(AvailCellsPerGroup / CellsPerGroup * 100, 2), "%"),
                    .groups = "drop"
                )
        })
    stats = bind_rows(stats)
}

# save the stats
write.table(
    stats,
    file = file.path(outdir, "stats.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
)

reporter$add(
    list(
        kind = "descr",
        content = "Statistics about the cells mapped to the comparison groups. Columns:"
    ),
    list(
        kind = "list",
        items = c(
            "_Group: The group name in the comparison, or null, if cells are not mapped to any group",
            "TotalCells: The total number of cells. This number should be the same for all groups",
            "CellsPerGroup: The number of cells in the mapped group",
            paste0(
                "AvailCellsPerGroup: The number of cells with CDR3 length between ",
                CDR3_MINLEN,
                " and ",
                CDR3_MAXLEN,
                " for each group. These cells are used for the analysis"
            ),
            "AvailCellsPct: The percentage of AvailCellsPerGroup over CellsPerGroup"
        )
    ),
    list(
        kind = "table",
        src = file.path(outdir, "stats.txt")
    ),
    h1 = "Available Cells"
)



log$info("Add amino acid features")
mdata = mdata %>%
    filter(!is.na(.Group) & length >= CDR3_MINLEN & length <= CDR3_MAXLEN) %>%
    add_percentAA() %>%
    add_positionalAA()


do_one_subset = function(s) {
    if (!is.null(s)) {
        log$info(paste("Processing subset", s))
    }
    if (is.null(s)) {
        data = mdata
        odir = file.path(outdir, "ALL")
    } else {
        data = mdata %>% filter(.Subset == s)
        odir = file.path(outdir, slugify(s))
    }
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    fits = list()
    for (len in CDR3_MINLEN:CDR3_MAXLEN){
        data_fit = data[data$length==len,]
        if (nrow(data_fit) == 0) { next }
        poss = get_positions(seq_length=as.numeric(as.character(len)))
        poss = paste("p", poss, sep="")
        feats = c("pI", "hydrophob", "volume")

        x <- matrix(nrow=nrow(data_fit), ncol=0)
        for (i in 1:length(poss)){
            ind = which(colnames(data_fit)==poss[i])
            for (j in 1:3){
            v = sapply(data_fit[,ind,drop=TRUE], function(z) get_feat_score(z, AA_MAPS[[j]]))
            vals = scale(v)[,1]
            vals[is.na(vals)] = 0
            x = cbind(x, vals)
            colnames(x)[ncol(x)] = paste(poss[i], feats[j], sep="_")
            }
        }
        y = ifelse(data_fit$.Group == target, 1, 0)
        if (any(table(y) <= 3) || length(table(y)) < 2) {
            if (is.null(s)) {
                log$warn(paste0("Not enough observations for target group '", target, "' with CDR3 length ", len, ". At least 4 observations are required."))
            } else {
                log$warn(paste0("Not enough observations for target group '", target, "' in subset '", s, "' with CDR3 length ", len, ". At least 4 observations are required."))
            }
        }
        # one multinomial or binomial class has 1 or 0 observations; not allowed
        if (any(table(y) <= 1)) { next }
        fit = glmnet(x, y, data=data_fit, alpha=0, lambda=0.01, family="binomial")
        fits[[len]] = tidy(fit)
    }


    # save the fits
    alldf = data.frame(matrix(nrow=0, ncol=4))
    colnames(alldf) = c("length", "imgt_pos", "feature", "estimate")
    for (len in CDR3_MINLEN:CDR3_MAXLEN){
        res = fits[[len]]
        if (is.null(res)) { next }
        res$length = len
        res$imgt_pos = sapply(res$term, function(x) strsplit(x, "_")[[1]][1])
        res$feature = sapply(res$term, function(x) strsplit(x, "_")[[1]][2])
        alldf = rbind(alldf, res[,colnames(alldf)])
    }
    alldf = alldf[alldf$feature %in% c("volume", "pI", "hydrophob"),]
    alldf$feature = as.character(alldf$feature)
    alldf$feature[alldf$feature=="hydrophob"] = "hydrophobicity"
    alldf$feature[alldf$feature=="pI"] = "isoelectric point"
    alldf$imgt_pos = as.character(alldf$imgt_pos)
    alldf$imgt_pos = factor(alldf$imgt_pos, levels=c("p104", "p105", "p106", "p107", "p108", "p109", "p110", "p111", "p111.1", "p112.1", "p112", "p113", "p114", "p115", "p116", "p117", "p118"))

    write.table(alldf, file = file.path(odir, "estimates.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    # save the plots
    gr <- alldf %>%
        group_by(imgt_pos, feature) %>%
        summarise(coef = mean(estimate))
    # Avoid too large values
    gr$coef[gr$coef > 1.5] <- 1.5
    gr$coef <- exp(gr$coef)  # Exponentiate the coefficients

    g <- LinePlot(gr, x = "imgt_pos", y = "coef", group_by = "feature",
        add_line = 1, x_text_angle = 90, xlab = "TCR position",
        ylab = paste("Coefficient for", target, "prediction"), title = s)

    save_plot(g, file.path(odir, "estimated_coefficients"),
        devpars = list(width = 1000, height = 1000, res = 100),
        formats = c("png", "pdf"))

    reporter$add(
        list(
            kind = "descr",
            content = "Estimated coefficients for each feature and position in the CDR3"
        ),
        h1 = ifelse(
            is.null(s),
            "Estimated OR (per s.d.)",
            paste0(paste(each_cols, collapse = ", "), " - ", s)
        ),
        h2 = ifelse(
            is.null(s),
            "#",
            "Estimated OR (per s.d.)"
        )
    )

    reporter$add(
        list(
            name = "Plot",
            contents = list(
                list(
                    kind = "image",
                    src = file.path(odir, "estimated_coefficients.png"),
                    download = file.path(odir, "estimated_coefficients.pdf")
                )
            )
        ),
        list(
            name = "Estimates",
            contents = list(
                list(
                    kind = "table",
                    src = file.path(odir, "estimates.txt")
                )
            )
        ),
        h1 = ifelse(
            is.null(s),
            "Estimated OR (per s.d.)",
            paste0(paste(each_cols, collapse = ", "), " - ", s)
        ),
        h2 = ifelse(
            is.null(s),
            "#",
            "Estimated OR (per s.d.)"
        ),
        ui = "tabs"
    )

    # distributions
    data$mid_hydro = sapply(data$midseq, function(x) get_feat_score(x, AA_MAPS[[2]]))
    data$smid_hydro = scale(data$mid_hydro)[,1]

    g <- RidgePlot(
        data = data,
        x = "smid_hydro",
        group_by = ".Group",
        xlab = "CDR3bmr hydrophobicity",
        ylab = "",
        add_vline = TRUE,
        alpha = 0.5,
        title = s,
        flip = TRUE
    )

    save_plot(g, file.path(odir, "distribution"),
        devpars = list(width = 1000, height = 1000, res = 100),
        formats = c("png", "pdf"))

    reporter$add(
        list(
            kind = "table_image",
            descr = paste0(
                "Hydrophobicity values are averaged over the CDR3 for each TCR and ",
                "then scaled to have a mean of 0 and a variance of 1. ",
                "Horizontal lines depict the mean for each population"
            ),
            src = file.path(odir, "distribution.png"),
            download = file.path(odir, "distribution.pdf")
        ),
        h1 = ifelse(
            is.null(s),
            "Hydrophobicity Distribution",
            paste0(paste(each_cols, collapse = ", "), " - ", s)
        ),
        h2 = ifelse(
            is.null(s),
            "#",
            "Hydrophobicity Distribution"
        )
    )

}

if (is.null(each_cols)) {
    do_one_subset(NULL)
} else {
    subsets = na.omit(unique(obj$.Subset))
    sapply(subsets, do_one_subset)
}

reporter$save(joboutdir)
