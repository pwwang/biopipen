library(methods)
library(fastLiquidAssociation)
{{rimport}}('__init__.r', 'plot.r')
options(stringsAsFactors = FALSE)
set.seed(8525)

infile   = {{i.infile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
casefile = {{i.casefile | R}}
inopts   = {{args.inopts | R}}
plotla   = {{args.plot | R}}
pval     = {{args.pval | R}}
dofdr    = {{args.fdr | R}}
fdrfor   = {{args.fdrfor | R}}
devpars  = {{args.devpars | R}}
nthread  = {{args.nthread | R}}
zcat     = {{args.zcat | R}}
ggs      = {{args.ggs | R}}

if (dofdr == T) dofdr = 'BH'

onecase = function(idata, zs, xs = NULL, nthread, zcat, case) {
	ncut   = max(ceiling(nrow(idata)/22), 4)
	inames = colnames(idata)
	z      = na.omit(match(zs, inames))
	x      = ifelse(is.null(xs), NULL, na.omit(match(xs, inames)))
	if (length(z) == 0)
		stop('Find no rows in "Z" group from input file.')
	ret    = fastMLA(
		as.matrix(idata), 
		nvec    = list(z = z, x = x),
		rvalue  = 0.3,
		cut     = ncut,
		zcat    = zcat,
		threads = nthread
	)
	if (nrow(ret) == 0) {
		return(cbind(ret, Pval = double()))
	} else {
		dupFlag = F
		if (nrow(ret) == 1) {
			dupFlag = T
			ret = rbind(ret, ret)
		}
		ret = mass.CNM(data=idata,GLA.mat=ret,nback=max(2, nrow(ret)))$`top p-values`
		ret = ret[, c(1:5, 9), drop = F]
		colnames(ret) = c('X', 'Y', 'Z', 'rhoDiff', 'LA', 'Pval')
		if (dupFlag)
			ret = ret[1,,drop=F]
		ret = cbind(Case = case, ret)
		return(ret)
	}
}

indata = t(read.table.inopts(infile, inopts))
ret = data.frame(
	Case    = character(),
	X       = character(),
	Y       = character(),
	Z       = character(),
	rhoDiff = double(),
	LA      = double(),
	Pval    = double()
)
if (dofdr != FALSE) {
	ret = cbind(ret, Qval = double())
}

casedata = read.table(casefile, header = F, row.names = NULL, check.names = F, sep = "\t")
if (ncol(casedata) == 2)
	casedata$Case = 'Case1'

cases = levels(factor(casedata[,3]))
for (case in cases) {
	cdata = casedata[which(casedata[,3] == case),,drop=F]
	zs    = cdata[which(cdata[,2] == 'Z'), 1]
	xs    = cdata[which(cdata[,2] == 'X'), 1]
	if (length(xs) == 0 || is.na(xs)) xs = NULL
	
	out = onecase(indata, zs, xs, nthread, zcat, case)
	

	if (dofdr != FALSE && fdrfor == 'case') {
		out = cbind(out, Qval = p.adjust(out$Pval, method = dofdr))
	}
	ret = rbind(ret, out)
}
if (dofdr != FALSE && fdrfor == 'all') {
	ret = cbind(ret, Qval = p.adjust(ret$Pval, method = dofdr))
}
ret = ret[ret$Pval < pval,,drop = FALSE]
write.table(pretty.numbers(ret, list(rhoDiff..LA = '%.3f', Pval..Qval = '%.2E')), outfile, col.names = T, row.names = F, sep = "\t", quote = F)

# plot
if (plotla && nrow(ret)>0) {
	ncut = max(ceiling(nrow(indata)/22), 4)
	for (i in 1:nrow(ret)) {
		row = ret[i,]
		plotfile = file.path(outdir, sprintf('%s.%s.', case, paste(row[2:4], collapse = '_')))
		if (!zcat) {
			LiquidAssociation::plotGLA(indata[, unlist(row[2:4]), drop = F], cut = ncut, dim = 3, pch = 4, filen = plotfile, save = T)
		} else {
			plotdata = data.frame(
				X = indata[, row$X], 
				Y = indata[, row$Y],
				Z = as.character(indata[, row$Z]))
			rcase = make.names(row$Z)
			colnames(plotdata) = c(row$X, row$Y, rcase)
			
			subgroups = levels(factor(plotdata[,3]))
			labels = sapply(subgroups, function(m) {
				r = cor(plotdata[plotdata[,3]==m, 1], plotdata[plotdata[,3]==m, 2], use = "pairwise.complete.obs")
				paste0(row$Z, '=', m, ' (r=', round(r, 3), ')')
			})
			ggs1 = c(ggs, list(
				geom_smooth = list(
					aes_string(color = rcase),
					method  = 'lm',
					se      = F
				),
				theme = list(
					legend.position = "bottom"
				),
				guides = list(
					color=guide_legend(ncol=1),
					shape=F
				),
				scale_color_manual = list(
					values = scales::hue_pal()(length(subgroups)), 
					name   = "",
					limit  = subgroups,
					labels = labels
				)
			))
			plot.scatter(
				plotdata, 
				paste0(plotfile, '.png'), 
				x      = 1,
				y      = 2,
				ggs    = ggs1,
				params = list(aes_string(shape = rcase, color = rcase))
			)
		}
	}
}



