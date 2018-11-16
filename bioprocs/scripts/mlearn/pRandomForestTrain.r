library(randomForest)
{{rimport}}('__init__.r')

infile       = {{i.infile | quote}}
outmodel     = {{o.outmodel | quote}}
outdir       = {{o.outdir | quote}}
args.plot    = {{args.plot | R}}
args.formula = {{args.formula | R}}
inopts       = {{args.inopts | R}}
devpars      = {{args.devpars | R}}
na           = {{args.na | R}}
prefix       = {{i.infile | stem | quote}}

indata = read.table.inopts(infile, inopts)
indata[is.na(indata)] = na

cnames = colnames(indata)
if (is.null(cnames)) {
	colnames(indata) = c(paste('X', 1:(ncol(indata)-1)), 'Y')
	cnames = colnames(indata)
}

if (is.null(args.formula)) {
	backtick = function(x) sprintf('`%s`', x)
	ycol = cnames[ncol(indata)]
	args.formula = as.formula(paste(cnames[ncol(indata)], '~', paste(backtick(cnames[1:(ncol(indata)-1)]), collapse = '+')))
} else {
	ycol = unlist(strsplit(args.formula, '\\s*~\\s*'))[1]
	args.formula = as.formula(args.formula)
}

model.rf = randomForest(args.formula, data = indata)
saveRDS(model.rf, outmodel)

# save importance
write.table(
	round(model.rf$importance[order(model.rf$importance[, 1], decreasing = T), , drop = F], 3), 
	file.path(outdir, paste0(prefix, '.importance.txt')), 
	col.names = T, row.names = T, quote = F, sep = "\t")

# plot importance
{% if args.plot %}
do.call(png, c(list(file.path(outdir, paste0(prefix, '.importance.png'))), devpars))
randomForest::varImpPlot(model.rf, main = 'Feature importance')
dev.off()
{% endif %}


