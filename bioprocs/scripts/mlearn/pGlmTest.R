{{rimport}}('__init__.r', 'plot.r')

infile  = {{i.infile | quote}}
inmodel = {{i.model | quote}}
outdir  = {{o.outdir | quote}}
outprob = {{args.outprob | R}}
inopts  = {{args.inopts | R}}
outauc  = {{args.outauc | R}}
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}
prefix  = file.path(outdir, {{i.infile | stem | quote}})

indata = read.table.inopts(infile, inopts)

cnames = colnames(indata)
if (is.null(cnames)) {
	colnames(indata) = c(paste('X', 1:(ncol(indata)-1)), 'Y')
	cnames = colnames(indata)
}

model = readRDS(inmodel)
prob  = do.call(predict.glm, c(list(model, indata, type="response"), params))

if (model$yval == 'prob') {
	if (outprob) {
		out = cbind(indata, Predict.Result = round(prob, 3), Predict.Prob = prob)
	} else {
		out = cbind(indata, Predict.Result = round(prob, 3))
	}
} else if (model$yval == 'categ') {
	ymin  = model$ymin
	ymax  = model$ymax
	index = round(prob * (ymax - ymin) + ymin)
	if (outprob) {
		out = cbind(indata, Predict.Result = model$ylabs[index], Predict.Prob = round(prob, 3))
	} else {
		out = cbind(indata, Predict.Result = model$ylabs[index])
	}
} else if (model$yval == 'numeric') {
	ymin  = model$ymin
	ymax  = model$ymax
	value = prob * (ymax - ymin) + ymin
	if (outprob) {
		out = cbind(indata, Predict.Result = value, Predict.Prob = round(prob, 3))
	} else {
		out = cbind(indata, Predict.Result = value)
	}
}
if (outprob) {
	preds = predict.glm(model, newdata = indata, type = "link", se.fit = TRUE)
	ci95_1 = model$family$linkinv(preds$fit - 1.96*preds$se.fit)
	ci95_2 = model$family$linkinv(preds$fit + 1.96*preds$se.fit)
	out = cbind(out, Predict.ProbCI95_1 = round(ci95_1, 3), Predict.ProbCI95_2 = round(ci95_2, 3))
}

ycol = model$ycol
if (outauc) {
	if (!ycol %in% cnames) {
		stop('No Y values found in input data, unable to plot ROC or output AUC.')
	}
	if (model$yval == 'prob' || model$yval == 'numeric') {
		stop('ROC/AUC is only available for categorical Y values, unable to plot ROC or output AUC.')
	}
	aucdata = data.frame(D = indata[, ycol] == out$Predict.Result, M = prob)
	roc = plot.roc(aucdata,
		paste0(prefix, '.roc.png'),
		params = list(returnTable = TRUE),
		ggs = ggs, devpars = devpars)
	out = cbind(out, AUC = roc$table[1, 'auc'])
}

write.table(out,
	paste0(prefix, '.result.txt'),
	row.names = list.get(inopts, 'rnames', TRUE),
	col.names = TRUE, sep = "\t", quote = F)
