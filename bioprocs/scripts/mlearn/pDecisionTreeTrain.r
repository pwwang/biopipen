library(rpart)
{{rimport}}('__init__.r', 'plot.r')

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

model.dt = rpart(args.formula, data = indata)
saveRDS(model.dt, outmodel)

# save importance
write.table(
	round(model.dt$variable.importance, 3), 
	file.path(outdir, paste0(prefix, '.importance.txt')), 
	col.names = F, row.names = T, quote = F, sep = "\t")

# plot importance
{% if args.plot %}
vi    = as.matrix(model.dt$variable.importance)
nfeat = min(30, nrow(vi))
vi    = vi[1:nfeat, , drop = F]
feats = rownames(vi)
plot.col(
	data.frame(Feature = feats, Importance = vi[, 1]), 
	file.path(outdir, paste0(prefix, '.importance.png')), 
	params  = list(aes(x = factor(feats, levels = rev(feats)))),
	ggs = list(coord_flip = list()),
	devpars = devpars)
{% endif %}


