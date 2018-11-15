{{rimport}}('__init__.r', 'plot.r')

infile       = {{i.infile | quote}}
outmodel     = {{o.outmodel | quote}}
outdir       = {{o.outdir | quote}}
args.plot    = {{args.plot | R}}
args.formula = {{args.formula | R}}
ggs          = {{args.ggs | R}}
devpars      = {{args.devpars | R}}
prefix       = {{i.infile | stem | quote}}
yval         = {{args.yval | quote}}

indata = read.table.nodup(
	infile, 
	header      = {{args.inopts.get('cnames', True) | R}},
	row.names   = {{args.inopts.get('rnames', True) | :1 if a else None | R}},
	sep         = {{args.inopts.get('delimit', '\t') | R}},
	check.names = {{args.inopts.get('check.names', False) | R}}
)

cnames = colnames(indata)
if (is.null(cnames)) {
	colnames(indata) = c(paste('X', 1:(ncol(indata)-1)), 'Y')
	cnames = colnames(indata)
}

if (is.null(args.formula)) {
	ycol = cnames[ncol(indata)]
	args.formula = as.formula(paste(cnames[ncol(indata)], '~', paste(cnames[1:(ncol(indata)-1)], collapse = '+')))
} else {
	ycol = unlist(strsplit(args.formula, '\\s*~\\s*'))[1]
	args.formula = as.formula(args.formula)
}

if (yval == 'categ') {
	ylabs   = levels(as.factor(indata[, ycol]))
	ylevels = order(ylabs)
	ymin    = min(ylevels)
	ymax    = max(ylevels)
	indata[, ycol] = (match(indata[, ycol], ylabs) - ymin)/(ymax - ymin)
	model.fit = lm(args.formula, data = indata)
	model.fit$ymin    = ymin
	model.fit$ymax    = ymax
	model.fit$ylabs   = ylabs
} else if (yval == 'numeric') {
	ymin = min(indata[, ycol])
	ymax = max(indata[, ycol])
	indata[, ycol] = (indata[, ycol] - ymin) / (ymax - ymin)
	model.fit = lm(args.formula, data = indata)
	model.fit$ymin = ymin
	model.fit$ymax = ymax
}
model.fit$yval = yval
model.fit$ycol = ycol

saveRDS(model.fit, outmodel)

model.fit.sum = summary(model.fit)
ret = model.fit.sum$coefficients

# confident interval
conf.int = confint(model.fit)
colnames(conf.int) = c('confInt25', 'confInt975')

OR = as.matrix(coef(model.fit))
colnames(OR) = 'OR'
ret = cbind(ret, exp(OR), exp(conf.int))
ret[, c(1:3, 5:7)] = round(ret[, c(1:3, 5:7)], 3)
ret[, 4] = formatC(ret[, 4], format = "E", digits = 2)
write.table(ret, file.path(outdir, paste0(prefix, '.features.txt')), row.names = T, col.names = T, sep = "\t", quote = F)

{% if args.plot %}
library(ggplot2)
plotfile = file.path(outdir, paste0(prefix, '.lm.png'))
plotdata = cbind(model.fit$model[, 1, drop=F], apply(model.fit$model[, -1, drop = F], 1, sum))
colnames(plotdata) = c(ycol, if (ncol(model.fit$model) > 2) 'Features' else colnames(model.fit$model)[ncol(model.fit$model)])

plot.scatter(plotdata, plotfile, x = 2, y = 1, ggs = ggs, devpars = devpars)
{% endif %}


