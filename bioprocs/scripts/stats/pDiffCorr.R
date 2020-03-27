
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
method   = {{args.method  | quote}}
ggs      = {{args.ggs     | R}}
devpars  = {{args.devpars | R}}
nthread  = {{args.nthread | R}}
stacked  = {{args.stacked | ?isinstance: dict | !:dict(row=_, col=_) | $R}}
if (dofdr == TRUE) dofdr = 'BH'
if (is.false(colfile)) {
    stop("A file must be specified to define groups in columns.")
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
            out = levels(as.factor(col))
            if (length(out) != 2) {
                stop(paste('Exact 2 groups needed for column group', name))
            }
            out
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
read_row = function(rowfile, stked) {
    if (is.false(rowfile)) {
        list(group=data.frame(RowGroup = c('Group1', 'Group2')),
             data=data.frame(
            ROWNAME=rep(inrows, 2),
            Group=rep(c('Group1', 'Group2'), each=length(inrows))
        ))
    } else {
        read_col(rowfile, stked)
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
                        Case = paste("Case", allgroups[i],
                                     allgroups[j], sep="_"),
                        SampleGroups = paste(allgroups[c(i, j)], collpase = ":")
                    )
                }
            }
        } else {
            for (group in names(coldata$group)) {
                casedata[nrow(casedata)+1, ] = list(
                    Case = paste("Case", group, sep="_"),
                    SampleGroups = paste(coldata$group[1:2, group],
                                         collapse = ':')
                )
            }
        }
    } else {
        casedata = read.table.inopts(casefile, list(rnames=FALSE, cnames=FALSE))
        casedata = row.apply(casedata, cnames = NULL, function(row) {
            row = unlist(row)
            if (!grepl(":", row[2], fixed = TRUE)) {
                if (is.null(coldata$group[[row[2]]])) {
                    stop(paste("Column group is not defined:", row[2]))
                }
                row[2] = paste(coldata$group[[row[2]]], collapse = ":")
            }
            row
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
                        one
                    } else {
                        paste(rowdata$group[[one]], collapse = ":")
                    }
                })
            }
        })
    }

    row.apply(ret, cnames = c("Case", "SampleGroup1", "SampleGroup2",
                              "RowGroup1", "RowGroup2"),
              function(row) {
                  row = unlist(row)
                  c(row[1],
                    unlist(strsplit(as.character(row[2]), ":", fixed = TRUE)),
                    unlist(strsplit(as.character(row[3]), ":", fixed = TRUE)))
              })
}

cases = read_case(casefile)

diff_corr = function(corr1, corr2, n1, n2) {
    # returns:
    #   Elem1  Elem2 Corr1 Corr2 N1 N2 Z       Pval
    #   A      B     .9    .3     10 20 1.435 .001
    corrs = cbind(melt(corr1), melt(corr2))[, c(1:3, 6), drop=FALSE]
    colnames(corrs) = c("Elem1", "Elem2", "Corr1", "Corr2")
    corrs$N1 = n1
    corrs$N2 = n2
    corrs$Z = apply(corrs, 1, function(row) {
        z1 = log( (1+as.numeric(row[3]))/(1-as.numeric(row[3])) ) / 2
        z2 = log( (1+as.numeric(row[4]))/(1-as.numeric(row[4])) ) / 2
        (z1 - z2) / sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
    })
    corrs$Pval = 2*pnorm(abs(corrs$Z), lower.tail = FALSE)
    corrs
}

get_corr = function(samples, rowgroup1, rowgroup2) {
    # get corr between elements in rowgroup1 and rowgroup2 using samples
    data1 = t(indata[rowdata$data[rowdata$data[, 2] == rowgroup1, 1], samples])
    data2 = t(indata[rowdata$data[rowdata$data[, 2] == rowgroup2, 1], samples])
    cor(data1, data2, method=method)
}

do_one_case = function(i) {
    case = cases[i, 1]
    samples1 = coldata$data[coldata$data[,2] == cases[i, 2], 1, drop = TRUE]
    samples2 = coldata$data[coldata$data[,2] == cases[i, 3], 1, drop = TRUE]
    rowgroup1 = cases[i, 4]
    rowgroup2 = cases[i, 5]

    corr1 = get_corr(samples1, rowgroup1, rowgroup2)
    corr2 = get_corr(samples2, rowgroup1, rowgroup2)

    cbind(data.frame(Case = case),
          diff_corr(corr1, corr2, length(samples1), length(samples2)))
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

if (is.true(doplot) && nrow(results)>0) {
    {{"plot.R" | rimport}}
    plot_one_case = function(i) {
        case = cases[i, 'Case']
        samgroup1 = cases[i, 2]
        samgroup2 = cases[i, 3]
        samples1 = coldata$data[coldata$data[,2] == cases[i, 2], 1, drop = TRUE]
        samples2 = coldata$data[coldata$data[,2] == cases[i, 3], 1, drop = TRUE]
        rowgroup1 = cases[i, 4]
        rowgroup2 = cases[i, 5]

        for (gene1 in rowdata$data[rowdata$data[, 2] == rowgroup1, 1]) {
            for (gene2 in rowdata$data[rowdata$data[, 2] == rowgroup2, 1]) {
                plotdata1 = as.data.frame(t(indata[c(gene1, gene2), samples1]))
                plotdata1$Group = samgroup1
                plotdata2 = as.data.frame(t(indata[c(gene1, gene2), samples2]))
                plotdata2$Group = samgroup2
                plotdata = rbind(plotdata1, plotdata2)
                rm(plotdata1, plotdata2)
                plotfile = file.path(outdir, paste0(case, "-", gene1,
                                                    "-", gene2, ".png"))
                result = results[results$Case == case &
                                 results$Elem1 == gene1 &
                                 results$Elem2 == gene2, , drop = FALSE]

                if (nrow(result) == 0 || result$Pval >= pcut) {next}

                plot.scatter(
                    plotdata, plotfile,
                    x = 1, y = 2,
                    ggs = c(ggs, list(
                        geom_smooth = list(
                            aes_string(color = "Group"),
                            method = 'lm',
                            se = FALSE
                        ),
                        scale_color_manual = list(
                            values = scales::hue_pal()(2),
                            limit  = c(samgroup1, samgroup2),
                            labels = c(
                                paste0(samgroup1,
                                       ' (r=', round(result$Corr1, 3),')'),
                                paste0(samgroup2,
                                       ' (r=', round(result$Corr2, 3),')')
                            )
                        ),
                        guides = list(shape = F)
                    )),
                    devpars = devpars,
                    params = list(aes_string(shape = "Group", color = "Group"))
                )
            }
        }
    }
    mcmapply(plot_one_case, 1:nrow(cases), mc.cores=nthread)
}
