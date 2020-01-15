{{'__init__.r', 'plot.r' | *rimport}}
library(randomForest)
library(dplyr)

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
indata = indata %>% mutate_if(is.character, as.factor)
for (col in cnames) {
	if (is.factor(model$call$data[, col])) {
		levels(indata[, col]) = levels(model$call$data[, col])
	}
}
prob  = do.call(predict, c(list(model, indata, type="response"), params))

if (outprob) {
	probs = do.call(predict, c(list(model, indata, type="prob")))
	out = cbind(indata, Predict.Result = as.vector(prob), Predict.Prob = probs[,1])
} else {
	out = cbind(indata, Predict.Result = as.vector(prob))
}

ycol = model$ycol
if (outauc) {
	if (!ycol %in% cnames) {
		stop('No Y values found in input data, unable to plot ROC or output AUC.')
	}
	if (model$type != 'classification') {
		stop('ROC/AUC is only available for categorical Y values, unable to plot ROC or output AUC.')
	}
	aucdata = data.frame(D = indata[, ycol], M = out$Predict.Prob)
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
