{{rimport}}('__init__.r', 'plot.r')
options(stringsAsFactors = FALSE)

infile  = {{i.infile | R}}
inmodel = {{i.model | R}}
outfile = {{o.outfile | R}}
outdir  = {{o.outdir | R}}
inopts  = {{args.inopts | R}}
outopts = {{args.out | R}}
plotroc = {{args.plot | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}

predictOne = function(idata, model, case) {
	prob  = predict(model, idata, type="response")
	ncols = ncol(idata)
	if (!is.null(model$ylabs)) { # categorical
		ymin   = model$ymin
		ymax   = model$ymax
		index  = round(prob * (ymax - ymin) + ymin)
		result = model$ylabs[index]
	} else {
		ymin   = model$ymin
		ymax   = model$ymax
		result = prob * (ymax - ymin) + ymin
	}
	out   = cbind(idata, Case = case, Predict.Result = result)
	colnames(out)[ncols+2] = paste0(model$ycol, '.Result')
	if (outopts$prob) {
		out = cbind(out, Predict.Prob = round(prob, 3))
		colnames(out)[ncols+3] = paste0(model$ycol, '.Prob')
	}
	out
}

indata = read.table.inopts(infile, inopts)
if (dir.exists(inmodel)) {
	mfiles = Sys.glob(file.path(inmodel, '*.*.rds'))
} else {
	mfiles = inmodel
}

out   = NULL
ycol  = NULL
ylabs = NULL
for (mfile in mfiles) {
	case  = tools::file_path_sans_ext(basename(mfile))
	model = readRDS(mfile)
	if (!is.null(ycol) && model$ycol != ycol)
		stop('Models are not on the same column')
	if (is.null(ycol)) ycol = model$ycol
	if (is.null(ylabs)) ylabs = model$ylabs
	ret   = predictOne(indata, model, case)
	out   = ifelse(is.null(out), ret, rbind(out, ret))
}

vars = colnames(indata)
if (!ycol %in% vars) {
	logger('Y column not included, no ROC plotted and no AUCs output anyway.')
	quit()
}

if (outopts$auc || is.list(plotroc)) {
	cases    = levels(factor(out$Case))
	if (length(ylabs) == 0) {
		rocdata  = data.frame(D = as.integer(out[, ycol] > mean(out[, ycol])))
	} else if (length(ylabs) == 1) {
		stop('Unable to do AUC on 1 level only for Y column.')
	} else if (length(ylabs) == 2) {
		rocdata = data.frame(D = as.integer(out[, ycol] == ylabs[2]))
	} else {
		rocdata = data.frame(D = as.integer(out[, ycol] == ylabs[length(ylabs)]))
	}
	for (case in cases) {
		rocdata  = cbind(rocdata, x = out[which(out$Case == case), paste0(ycol, '.Result')])
		colnames(rocdata)[ncol(rocdata)] = case
	}
	nrows = nrow(rocdata)
	aucs = plot.roc(
		rocdata, 
		plotfile = file.path(outdir, {{i.infile | fn2 | @append: '.roc.png' | quote}}),
		stacked  = FALSE,
		params   = c(list(returnAUC = TRUE), plotroc),
		ggs      = ggs,
		devpars  = devpars
	)
	if (outopts$auc) {
		out$AUC = NA
		for (case in cases)
			out[which(out$Case == case), 'AUC'] = aucs[[case]]
	}
}

write.table(out, outfile, row.names = inopts$rnames, col.names = T, sep = "\t", quote = F)
