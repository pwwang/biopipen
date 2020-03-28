
{{"__init__.R" | rimport}}

library(methods)
library(reshape2)
library(dplyr)
library(parallel)
options(stringsAsFactors = FALSE)

infile   = {{i.infile     | quote}}
colfile  = {{i.colfile    | quote}}
casefile = {{i.casefile   | quote}}
rowfile  = {{i.rowfile    | quote}}
outdir   = {{o.outdir     | quote}}
outfile  = {{o.outfile    | quote}}
inopts   = {{args.inopts  | R}}
pcut     = {{args.pval    | R}}
dofdr    = {{args.fdr     | R}}
doplot   = {{args.plot    | R}}
ggs      = {{args.ggs     | R}}
devpars  = {{args.devpars | R}}
nthread  = {{args.nthread | R}}
stacked  = {{args.stacked | ?isinstance: dict | !:dict(row=_, col=_) | $R}}
if (dofdr == TRUE) dofdr = 'BH'
if (is.false(colfile)) {
    stop("A file must be specified to define groups in columns.")
}
if (doplot) {
    {{"plot.R" | rimport}}
}

# read colfile and rowfile into:
# list(
#	group =
#   	+------------+------------+------------+------+-----------+
#       | SamGroup1  | SamGroup2  | SamGroup3  | ... |	SamGroupN |
#   	+------------+------------+------------+-----+------------+
#       | SubGroup11 | SubGroup21 | SubGroup31 | ... | SubGroupN1 |
#       | SubGroup12 | SubGroup22 | SubGroup32 | ... | SubGroupN2 |
#   	+------------+------------+------------+-----+------------+
#   data =
#       +---------+------------+
#       | Sample  | Group      |
#       +---------+------------+
#       | Sample1 | SubGroup11 |
#       | Sample1 | SubGroup21 |
#       | Sample2 | SubGroup22 |
#       | ...     | ...        |
#       | SampleX | SubGroupN2 |
#       +---------+------------+
# )
read_col = function(colfile, stked) {
    if (stked) {
        coldata = read.table.inopts(colfile, list(cnames=FALSE, rnames=FALSE))
        list(group = NULL, data = coldata)
    } else {
        coldata = read.table.inopts(colfile, list(cnames=TRUE, rnames=TRUE))
        # attach the column name to the groups
        coldata = col.apply(coldata, func = function(col, name) {
            col = as.vector(col)
            col[!is.na(col)] = paste0(name, '-', na.omit(unlist(col)))
            factor(col, levels = unique(col))
        })
        retgroup = col.apply(coldata, rnames = NULL, function(col, name) {
            out = na.omit(unique(unlist(col)))
            lout = length(out)
            if (lout < 2) {
                stop(paste('At least 2 groups needed for column group', name))
            }
            c(out, rep(NA, 10-lout))
        })

        cnames = paste0("V", 1:ncol(coldata))
        colnames(coldata) = cnames
        coldata$ROWNAME = rownames(coldata)
        retdata = melt(data = coldata, id.vars = "ROWNAME",
                       measure.vars = cnames, na.rm = TRUE)
        list(group = as.data.frame(retgroup),
             data = retdata %>% select(ROWNAME, value))
    }
}

indata = read.table.inopts(infile, inopts)
inrows = rownames(indata)
incols = colnames(indata)
indata = as.data.frame(t(indata))
read_row = function(rowfile, stked) {
    if (is.false(rowfile)) {
        list(group=data.frame(RowGroup = c('Group1', 'Group2')),
             data=data.frame(
            ROWNAME=rep(inrows, 2),
            Group=rep(c('Group1', 'Group2'), each=length(inrows))
        ))
    } else {
        ret = read_col(rowfile, stked)
        if (!stked) {
            ret$group = ret$group[1:2, , drop=FALSE]
        }
        ret
    }
}
coldata = read_col(colfile, stacked$col)
rowdata = read_row(rowfile, stacked$row)

read_case = function(casefile) {
    if (is.false(casefile)) {
        casedata = data.frame(Case = character(),
                              SampleGroups = character(),
                              stringsAsFactors = FALSE)
        if (is.null(coldata$group)) {
            allgroups = levels(coldata$data[,1])
            for (i in 1:(length(allgroups)-1)) {
                for (j in (i+1):length(allgroups)) {
                    casedata[nrow(casedata)+1, ] = list(
                        Case = paste(allgroups[i], allgroups[j], sep=":"),
                        SampleGroups = paste(allgroups[c(i, j)], collpase = ":")
                    )
                }
            }
        } else {
            for (group in names(coldata$group)) {
                casedata[nrow(casedata)+1, ] = list(
                    Case = group,
                    SampleGroups = paste(coldata$group[1:2, group],
                                         collapse = ':')
                )
            }
        }
    } else {
        casedata = read.table.inopts(casefile, list(rnames=FALSE, cnames=FALSE))
        casedata = row.apply(casedata, cnames = NULL, function(row) {
            row = unlist(row)
            case = row[1]
            if (!grepl(":", case, fixed = TRUE)) {
                if (is.null(coldata$group[[case]])) {
                    stop(paste("Column group is not defined:", case))
                }
                row[1] = paste(na.omit(coldata$group[[case]]), collapse = ":")
            }
            c(case, row)
        })
    }

    # no row groups specified
    if (ncol(casedata) == 2) {
        # exhaust all row (group) pairs
        casedata$RowGroups = ""
        if (is.null(rowdata$group)) {
            # exhaust all combinations
            groups = na.omit(unique(rowdata$data[,2]))
            ngroups = length(groups)
            ngpairs = (ngroups * (ngroups - 1)) / 2
            ret = col.apply(casedata, rnames = NULL, function(col, name) {
                if (name != "RowGroups") {
                    rep(col, each = ngpairs)
                } else {
                    out = c()
                    for (i in 1:(ngroups-1)) {
                        for (j in (i+1):ngroups) {
                            out = c(out, paste(groups[c(i, j)], collapse = ":"))
                        }
                    }
                    rep(out, length(col))
                }
            })
        } else {
            # row groups defined, exhaust them
            groups = colnames(rowdata$group)
            ngroups = length(groups)
            ret = col.apply(casedata, rnames = NULL, function(col, name) {
                if (name != "RowGroups") {
                    rep(col, each = ngroups)
                } else {
                    out = c()
                    for (group in groups) {
                        out = c(out, paste(rowdata$group[, group],
                                collapse = ":"))
                    }
                    rep(out, length(col))
                }
            })
        }
    # row groups specified
    } else {
        ret = col.apply(casedata, rnames = NULL, function(col, index) {
            if (index < 3) {
                as.vector(col)
            } else {
                sapply(col, function(one) {
                    if (is.null(rowdata$group[[one]])) {
                        if (grepl(":", one, fixed = TRUE)) {
                            one
                        } else {
                            warning(paste("No such column group:", one))
                            NA
                        }
                    } else {
                        paste(rowdata$group[[one]], collapse = ":")
                    }
                })
            }
        })
    }

    row.apply(ret[complete.cases(ret),,drop=FALSE],
              cnames = c("Case", "SampleGroups",
                         "RowGroup1", "RowGroup2"),
              function(row) {
                  row = unlist(row)
                  rowgroups = unlist(strsplit(as.character(row[3]),
                                              ":", fixed = TRUE))
                  c(row[1], row[2], rowgroups)
              })
}

# "Case", "SampleGroups", "RowGroup1", "RowGroup2"
cases = read_case(casefile)

chow_test_one = function(groups, row1, row2) {
    fmula = as.formula(paste(row2, "~", row1))
    pooled_lm = tryCatch({lm(fmula, data = indata)}, error = function(e) NULL)
    group_lms = sapply(groups, function(grup) {
        list(lm(fmula, data = indata[coldata$data[ coldata$data[,2] == grup, 1],
                                     c(row1, row2), drop=FALSE ]))
    })
    pooled.ssr <- ifelse(is.null(pooled_lm), NA, sum(pooled_lm$residuals^2))
    subssr <- ifelse(is.false(group_lms, "any"), NA,
                     sum(sapply(group_lms, function(x) sum(x$residuals^2))))
    ngroups <- length(groups)
    K <- 2 #+ length(covs)
    J <- (ngroups - 1) * K
    DF <- nrow(indata) - ngroups * K
    FS <- (pooled.ssr - subssr) * DF / subssr / J
    list(pooled_lm = pooled_lm,
         group_lms = group_lms,
         Fstat = FS,
         Pval = pf(FS, J, DF))
}

total = nrow(indata)

plot_chow = function(pooled_lm, group_lms, case, row1, row2) {
    plotdata <- do.call(
        rbind,
        lapply(names(group_lms),
               function(m) data.frame(
                   group_lms[[m]]$model[, c(row1, row2), drop = FALSE],
                   group = m
               ))
    )
    colnames(plotdata)[3] <- case
    if (!is.null(ggs$scale_color_discrete)) {
        ggs$scale_color_discrete$name <- ifelse(
            is.function(ggs$scale_color_discrete$name),
            ggs$scale_color_discrete$name(case),
            case
        )
        ggs$scale_color_discrete$labels <- sapply(
            names(group_lms),
            function(x) {
                coeff <- as.list(group_lms[[x]]$coefficients)
                bquote(.(x):beta[.(row1)] == .(round(coeff[[row1]], 3)) ~
                       "," ~ epsilon == .(round(coeff[["(Intercept)"]], 3)))
            }
        )
    }

    if (!is.null(ggs$scale_shape_discrete)) {
        ggs$scale_shape_discrete$name <- ifelse(
            is.function(ggs$scale_shape_discrete$name),
            ggs$scale_shape_discrete$name(case),
            case
        )
        ggs$scale_shape_discrete$labels <- sapply(
            names(group_lms),
            function(x) {
                coeff <- as.list(group_lms[[x]]$coefficients)
                bquote(.(x):beta[.(row1)] == .(round(coeff[[row1]], 3)) ~
                       "," ~ epsilon == .(round(coeff[["(Intercept)"]], 3)))
        })
    }

    if (is.null(ggs$scale_color_discrete) && is.null(ggs$scale_shape_discrete)){
        ggs$scale_color_discrete <- list(
            name = case,
            labels = sapply(names(group_lms), function(x) {
                coeff <- as.list(group_lms[[x]]$coefficients)
                bquote(.(x):beta[.(row1)] == .(round(coeff[[row1]], 3)) ~
                       "," ~ epsilon == .(round(coeff[["(Intercept)"]], 3)))
            })
        )
    }

    if (is.null(ggs$scale_color_discrete))
        ggs$scale_color_discrete <- ggs$scale_shape_discrete
    if (is.null(ggs$scale_shape_discrete))
        ggs$scale_shape_discrete <- ggs$scale_color_discrete

    plot.points(
        plotdata,
        file.path(outdir, paste0(case, "-", row1, "-", row2, ".png")),
        x = 1, y = 2,
        params = list(aes_string(color = case, shape = case)),
        ggs = c(ggs, list(
            geom_smooth = list(aes_string(color = case),
                               method = "lm", se = FALSE)
        ))
    )
}

do_one_case = function(i) {
    case = cases[i, 1]
    samgroups = unlist(strsplit(cases[i, 2], ":", fixed=TRUE))
    rowgroup1 = cases[i, 3]
    rowgroup2 = cases[i, 4]
    rowgroup1 = rowdata$data[ rowdata$data[,2] == rowgroup1, 1 ]
    rowgroup2 = rowdata$data[ rowdata$data[,2] == rowgroup2, 1 ]

    ret = data.frame(Case=character(), Elem1=character(), Elem2=character(),
                     TotalN=numeric(), Ns=character(), Fstat=numeric(),
                     Pval=numeric())
    for (row1 in rowgroup1) {
        for (row2 in rowgroup2) {
            out = chow_test_one(samgroups, row1, row2)
            ret[nrow(ret)+1, ] = list(Case=case, Elem1=row1, Elem2=row2,
                                      TotalN=total, Ns=paste(sapply(
                                          samgroups,
                                          function(grup) {
                                              sum(coldata$data[,2] == grup)
                                          }
                                      ), collapse=":"),
                                      Fstat=out$Fstat, Pval=out$Pval)
            if (!is.na(out$Pval) && out$Pval < pcut) {
                plot_chow(out$pool_lm, out$group_lms, case, row1, row2)
            }
        }
    }
    ret
}

results = do.call(rbind, mcmapply(do_one_case, 1:nrow(cases),
                                  SIMPLIFY=FALSE,
                                  mc.cores=nthread))

if (is.true(dofdr) && nrow(results) > 0) {
    results$Qval = p.adjust(results$Pval, method = dofdr)
}
results = results[!is.na(results$Pval) & results$Pval < pcut, , drop = FALSE]
write.table(results, outfile, row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")
