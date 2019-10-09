{{rimport}}('__init__.r')

infile       = {{i.infile | quote}}
outmodel     = {{o.outmodel | quote}}
outdir       = {{o.outdir | quote}}
args.plot    = {{args.plot | R}}
args.formula = {{args.formula | R}}
inopts       = {{args.inopts | R}}
params       = {{args.params | R}}
ggs          = {{args.ggs | R}}
devpars      = {{args.devpars | R}}
prefix       = {{i.infile | stem | quote}}
yval         = {{args.yval | quote}}

indata = read.table.inopts(infile, inopts)

cnames = colnames(indata)
if (is.null(cnames)) {
	colnames(indata) = c(paste('X', 1:(ncol(indata)-1)), 'Y')
	cnames = colnames(indata)
}

if (is.null(args.formula)) {
	ycol = cnames[ncol(indata)]
	args.formula = as.formula(paste(
		cnames[ncol(indata)], '~',
		paste(bQuote(cnames[1:(ncol(indata)-1)]),
		collapse = '+')))
} else {
	args.formula = as.formula(args.formula)
	ycol = all.vars(args.formula)[1]
}

if (yval == 'categ') {
	ylabs   = levels(as.factor(indata[, ycol]))
	ylevels = order(ylabs)
	ymin    = min(ylevels)
	ymax    = max(ylevels)
	indata[, ycol] = (match(indata[, ycol], ylabs) - ymin)/(ymax - ymin)
	model.fit = do.call(glm, c(list(args.formula, data = indata), params))
	model.fit$ylabs   = ylabs
} else if (yval == 'numeric') {
	ymin = min(indata[, ycol])
	ymax = max(indata[, ycol])
	indata[, ycol] = (indata[, ycol] - ymin) / (ymax - ymin)
	model.fit = do.call(glm, c(list(args.formula, data = indata), params))
}
model.fit$ymin = ymin
model.fit$ymax = ymax
model.fit$yval = yval
model.fit$ycol = ycol

saveRDS(model.fit, outmodel)

model.fit.sum = summary(model.fit)
ret = model.fit.sum$coefficients

OR = as.matrix(coef(model.fit))
colnames(OR) = 'OR'
ret = cbind(ret, exp(OR))
ret[, 1:3] = round(ret[, 1:3], 3)
ret[, 4] = formatC(ret[, 4], format = "E", digits = 2)

tryCatch({
	# confident interval
	conf.int = confint(model.fit)
	colnames(conf.int) = c('ci_95_1', 'ci_95_2')
	ret = cbind(ret, exp(cont.int))
}, error = function(e) {})

rownames(ret) = nobQuote(rownames(ret))
write.table(ret,
	file.path(outdir, paste0(prefix, '.features.txt')),
	row.names = TRUE,
	col.names = TRUE,
	sep = "\t",
	quote = FALSE)

{% if args.plot %}
{{rimport}}('plot.r')
plotfile = file.path(outdir, paste0(prefix, '.glm.png'))
featdata = model.fit$model[, -1, drop = F]
# what about categorical features
numcols = unlist(lapply(featdata, is.numeric))

plotdata = cbind(
	model.fit$model[, 1, drop=F],
	apply(featdata[, numcols, drop = FALSE], 1, sum))
colnames(plotdata) = c(
	ycol,
	ifelse(ncol(model.fit$model) > 2,
		'Features',
		colnames(model.fit$model)[ncol(model.fit$model)]))

plot.scatter(plotdata, plotfile, x = 2, y = 1, ggs = ggs, devpars = devpars)
{% endif %}
