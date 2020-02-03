
{{"__init__.R", "plot.R", "sampleinfo.R" | rimport}}

{% from pyppl.utils import always_list %}
infile   <- {{i.infile | quote}}
gfile    <- {{i.gfile | quote}}
prefix   <- {{i.infile | fn2 | quote}}
outdir   <- {{o.outdir | quote}}
inopts   <- {{args.inopts | R}}
tsform   <- {{args.tsform | ?!:'NULL'}}
filter   <- {{args.filter | ?!:'NULL'}}
annocols <- {{args.annocols | ?=always_list | $R}}
plots    <- {{args.plot | R}}
ggs      <- {{args.ggs | R}}
params   <- {{args.params | R}}
devpars  <- {{args.devpars | R}}

expr = read.table.inopts(infile, inopts)
if (is.true(annocols)) {
    if (is.numeric(annocols)) {
        annocols = colnames(expr)[annocols]
    }
    expr = expr[, setdiff(colnames(expr), annocols), drop=FALSE]
}
if (is.function(filter)) {
    expr = filter(expr)
    outfile = file.path(outdir, basename(infile))
    write.table(expr, outfile, row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE)
}
if (is.function(tsform)) {
    expr = tsform(expr)
}

saminfo = NULL
groups = NULL
batches = NULL
if (gfile != "") {
    saminfo = SampleInfo2$new(gfile)
    groups  = saminfo$all.groups()
    batches = saminfo$all.batches()
}

if (is.true(plots$pca)) {
    library(PCAtools)
    pcafile = file.path(outdir, paste0(prefix, '.pca.png'))
    metadata = NULL
    samples = colnames(expr)
    if (!is.null(groups)) {
        tmp = sapply(samples,
                     function(sam) saminfo$sample.info(sam, 'Group'))
        if (is.null(metadata)) {
            metadata = data.frame(Group = tmp)
            rownames(metadata) = samples
        } else {
            metadata$Group = tmp
        }
    }
    if (!is.null(batches)) {
        tmp = sapply(samples,
                     function(sam) saminfo$sample.info(sam, 'Batch'))
        if (is.null(metadata)) {
            metadata = data.frame(Batch = tmp)
            rownames(metadata) = samples
        } else {
            metadata$Batch = tmp
        }
    }
    p = pca(expr, metadata = metadata, removeVar = 0.1)
    biplot_params = list(pcaobj=p)
    if (!is.null(groups)) {
        biplot_params$shape = 'Group'
        biplot_params$legendPosition = 'right'
    }
    if (!is.null(batches)) {
        biplot_params$colby = 'Batch'
        biplot_params$legendPosition = 'right'
    }
    p = do.call(biplot, biplot_params)
    save.plot(p, pcafile, devpars = devpars)
}

if (is.true(plots$boxplot)) {
    bpfile = file.path(outdir, paste0(prefix, '.boxplot.png'))
    plot.boxplot(expr, bpfile, stacked = F, devpars = devpars,
                 params = params$boxplot, ggs = ggs$boxplot)

    if (!is.null(groups)) {
        bpfile = file.path(outdir,paste0(prefix, '.group.boxplot.png'))
        bpdata = NULL
        for (group in groups) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Group', group),
                                         drop = TRUE]),
                Group = group
            )
            if (is.null(bpdata)) {
                bpdata = tmp
            } else {
                bpdata = rbind(bpdata, tmp)
            }
        }
        plot.boxplot(bpdata, bpfile, stacked = TRUE, devpars = devpars,
                     params = params$boxplot, ggs = ggs$boxplot)
    }

    if (!is.null(batches)) {
        bpfile = file.path(outdir,paste0(prefix, '.batch.boxplot.png'))
        bpdata = NULL
        for (batch in batches) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Batch', batch),
                                         drop = TRUE]),
                Batch = batch
            )
            if (is.null(bpdata)) {
                bpdata = tmp
            } else {
                bpdata = rbind(bpdata, tmp)
            }
        }
        plot.boxplot(bpdata, bpfile, stacked = TRUE, devpars = devpars,
                     params = params$boxplot, ggs = ggs$boxplot)
    }
}

if (is.true(plots$violin)) {
    vlfile = file.path(outdir, paste0(prefix, '.violin.png'))
    plot.violin(expr, vlfile, stacked = F, devpars = devpars,
                params = params$violin, ggs = ggs$violin)

    if (!is.null(groups)) {
        vlfile = file.path(outdir,paste0(prefix, '.group.violin.png'))
        vldata = NULL
        for (group in groups) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Group', group),
                                         drop = TRUE]),
                Group = group
            )
            if (is.null(vldata)) {
                vldata = tmp
            } else {
                vldata = rbind(vldata, tmp)
            }
        }
        plot.violin(vldata, vlfile, stacked = TRUE, devpars = devpars,
                    params = params$violin, ggs = ggs$violin)
    }

    if (!is.null(batches)) {
        vlfile = file.path(outdir,paste0(prefix, '.batch.violin.png'))
        vldata = NULL
        for (batch in batches) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Batch', batch),
                                         drop = TRUE]),
                Batch = batch
            )
            if (is.null(vldata)) {
                vldata = tmp
            } else {
                vldata = rbind(vldata, tmp)
            }
        }
        plot.violin(vldata, vlfile, stacked = TRUE, devpars = devpars,
                    params = params$violin, ggs = ggs$violin)
    }
}

if (is.true(plots$histogram)) {
    histfile = file.path(outdir, paste0(prefix, ".histo.png"))
    plot.histo(stack(as.data.frame(expr)), histfile, devpars = devpars,
               params = params$histogram, ggs = ggs$histogram)
    if (!is.null(groups)) {
        histfile = file.path(outdir, paste0(prefix, ".group.histo.png"))
        hisdata = NULL
        for (group in groups) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Group', group),
                                         drop = TRUE]),
                Group = group
            )
            if (is.null(hisdata)) {
                hisdata = tmp
            } else {
                hisdata = rbind(hisdata, tmp)
            }
        }
        plot.histo(hisdata,
                   histfile,
                   devpars = devpars,
                   params = c(list(aes_string(fill="Group"), alpha=.7),
                              params$histogram),
                   ggs = ggs$histogram)
    }

    if (!is.null(batches)) {
        histfile = file.path(outdir, paste0(prefix, ".batch.histo.png"))
        hisdata = NULL
        for (batch in batches) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Batch', batch),
                                         drop = TRUE]),
                Batch = batch
            )
            if (is.null(hisdata)) {
                hisdata = tmp
            } else {
                hisdata = rbind(hisdata, tmp)
            }
        }
        plot.histo(hisdata,
                   histfile,
                   devpars = devpars,
                   params = c(list(aes_string(fill="Batch"), alpha=.7),
                              params$histogram),
                   ggs = ggs$histogram)
    }
}

if (is.true(plots$qqplot)) {
    if (!is.null(groups)) {
        if (length(groups) != 2) {
            stop('Cannot plot QQ plot for more than 2 groups.')
        }
        qqfile = file.path(outdir, paste0(prefix, '.group.qq.png'))
        qqdata = NULL
        for (group in groups) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Group', group),
                                         drop=TRUE]),
                Group = group
            )
            if (is.null(tmp)) {
                qqdata = tmp
            } else {
                qqdata = rbind(qqdata, tmp)
            }
        }
        plot.qq(qqdata,
                qqfile,
                devpars = devpars,
                params = params$qqplot,
                ggs = c(list(labs=list(x=groups[2], y=groups[1])), ggs$qqplot))
    }

    if (!is.null(batches)) {
        if (length(batches) != 2) {
            stop('Cannot plot QQ plot for more than 2 batches.')
        }
        qqfile = file.path(outdir, paste0(prefix, '.batch.qq.png'))
        qqdata = NULL
        for (batch in batches) {
            tmp = data.frame(
                Expression = unlist(expr[, saminfo$get.samples('Batch', batch),
                                         drop=TRUE]),
                Batch = batch
            )
            if (is.null(tmp)) {
                qqdata = tmp
            } else {
                qqdata = rbind(qqdata, tmp)
            }
        }
        plot.qq(qqdata,
                qqfile,
                devpars = devpars,
                params = params$qqplot,
                ggs = c(list(labs=list(x=groups[2], y=groups[1])), ggs$qqplot))
    }
}
