library(methods)
library(data.table)
{{rimport}}('__init__.r', 'plot.r')
options(stringsAsFactors = FALSE)

infile   = {{i.infile     | quote}}
samfile  = {{i.samfile    | quote}}
casefile = {{i.casefile   | quote}}
gfile    = {{i.groupfile  | quote}}
outdir   = {{o.outdir     | quote}}
outfile  = {{o.outfile    | quote}}
inopts   = {{args.inopts  | R}}
pcut     = {{args.pval    | R}}
dofdr    = {{args.fdr     | R}}
fdrfor   = {{args.fdrfor  | R}}
doplot   = {{args.plot    | R}}
method   = {{args.method  | quote}}
ggs      = {{args.ggs     | R}}
devpars  = {{args.devpars | R}}
if (dofdr == T) dofdr = 'BH'

getCorr = function(gname1, mat1, gname2, mat2) {
	# returns:
	# $corr
	# 	Group1	Group2	Corr
	# 	A		B		.9
	# 	A		C		.1
	# $n
	# 	10
	rnames1 = rownames(mat1)
	rnames2 = rownames(mat2)
	mat     = merge(mat1, mat2, by = c())
	halfcol = ncol(mat) / 2
	rm(mat1, mat2)
	corr = data.frame(corr = apply(mat, 1, function(row) cor(
		row[1:halfcol],
		row[(halfcol+1):(2*halfcol)],
		method = method,
		use = "pairwise.complete.obs"
	)))
	corr = cbind(merge(rnames1, rnames2, by = c()), corr)
	colnames(corr) = c(gname1, gname2, 'Corr')
	
	lvl1 = levels(corr[,1])
	lvl2 = levels(corr[,2])
	levels(corr[,1]) = union(lvl1, lvl2)
	levels(corr[,2]) = union(lvl1, lvl2)
	list(corr = corr[which(corr[,1] != corr[,2]),, drop = F], n = halfcol)
}

diffCorr = function(corr1, corr2) {
	# returns:
	#   Group1 Group2 Corr1 Corr2 N1 N2 Z     Pval  Qval
	#   A      B     .9    .3     10 20 1.435 .001  .01
	n1    = corr1$n
	n2    = corr2$n
	colnames(corr1$corr)[3] = 'Corr1'
	colnames(corr2$corr)[3] = 'Corr2'
	corr  = cbind(corr1$corr$Corr1, corr2$corr$Corr2)
	dzcor = apply(corr, 1, function(row) {
		z1 = log( (1+row[1])/(1-row[1]) ) / 2
		z2 = log( (1+row[2])/(1-row[2]) ) / 2
		dz = z1 - z2
		sedz = sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
		dz/sedz
		# two-sized
		#list(Z = z, Pval = 2 * pnorm(abs(z), lower.tail=F))
	})
	pval = 2*pnorm(abs(dzcor), lower.tail = F)
	ret = cbind(corr1$corr, Corr2 = corr2$corr$Corr2, N1 = n1, N2 = n2, Z = dzcor, Pval = pval)
}

indata  = read.table.inopts(infile, inopts)
#    S1   S2   S3 ...
# A  1    2    1
# B  3    9    8
# ...
samdata = read.table(samfile, header = F, row.names = NULL, check.names = F)
# S1     G1
# S2     G1
# S3     G2
# ...
# Sm     G3
# Sn     G3
samgrps = factor(samdata[,2])
samgrps = relevel(samgrps, as.character(samdata[1,2]))
samgrps = levels(samgrps)
# G1, G2, G3

if (gfile != '') {
	rgdata  = read.table(gfile, header = F, row.names = 1, check.names = F)
	# A     TF
	# B     TF
	# C     TF
	# ...
	# X     Gene
	# Y     Gene
	# Z     Gene
	rgroups = factor(rgdata[,1])
	rgroups = levels(rgroups)
	if (length(rgroups)!=2) {
		stop('Group file should only contain 2 groups!')
	}
	rgname1 = rgroups[1] # TF
	rgname2 = rgroups[2] # Gene
	rgroup1 = rownames(rgdata[which(rgdata[,1] == rgname1),, drop = F])
	# A, B, C
	rgroup2 = rownames(rgdata[which(rgdata[,1] == rgname2),, drop = F])
	# X, Y, Z
} else {
	rgname1 = 'R1'
	rgname2 = 'R2'
	rgroup1 = rownames(indata)
	# A, B, C, ..., X, Y, Z
	rgroup2 = rgroup1
	# A, B, C, ..., X, Y, Z
}

corrs = list()
for (samgrp in samgrps) {
	sams = samdata[which(samdata[,2] == samgrp), 1]
	if (length(sams) < 3) next # not able to do correlation
	corrs[[samgrp]] = getCorr(rgname1, indata[rgroup1, sams, drop=F], rgname2, indata[rgroup2, sams,drop=F])
	# $corr
	# 	G1	G2	Corr
	# 	A	Z	.9
	# 	A	Z	.1
	# 	...
	# $n
	# 	10
}

if (casefile != '') {
	cases = read.table(casefile, header = F, row.names = NULL, sep = "\t", check.names = F)
	colnames(cases) = c('x', 'y')
	# x     y
	# G1    G3
	# G2    G3
} else {
	cases = merge(samgrps, samgrps, by = c())
	lvlx  = levels(cases$x)
	lvly  = levels(cases$y)
	levels(cases$x) = union(lvlx, lvly)
	levels(cases$y) = union(lvlx, lvly)
	cases = cases[which(cases$x != cases$y),, drop = F]
	# x     y
	# G1    G3
	# G2    G3
}

results = NULL
# returns:
# 	Case1   Case2  TF Gene Corr1 Corr2 N1 N2 Z     Pval  Qval
# 	G1      G3     A  B    .9    .3    10 20 1.435 .001  .01
for (i in 1:nrow(cases)) {
	x = as.vector(cases[i, 1])
	y = as.vector(cases[i, 2])
	if (is.null(corrs[[x]]) || is.null(corrs[[y]])) 
		next
	dc = diffCorr(corrs[[x]], corrs[[y]])
	if (nrow(dc) == 0)
		dc = cbind(Case1 = NULL, Case2 = NULL, dc)
	else
		dc = cbind(Case1 = x, Case2 = y, dc)
	results = if(is.null(results)) dc else rbind(results, dc)
}
if (dofdr != F && !is.null(results)) {
	results$Qval = NA
	if (fdrfor == 'all') {
		results$Qval = p.adjust(results$Pval, method = dofdr)
	} else {
		for (i in 1:nrow(cases)) {
			x = as.vector(cases[i, 1])
			y = as.vector(cases[i, 2])
			rows = which(results$Case1 == x & results$Case2 == y)
			if (length(rows)>0)
				results[rows, 'Qval'] = p.adjust(results[rows, 'Pval'], method = dofdr)
		}
	}
}
if (is.null(results)) {
	if (dofdr == F) {
		results = data.frame(Case1 = character(), Case2 = character(), rgname1 = character(), rgname2 = character(), Corr1 = double(), Corr2 = double(), N1 = integer(), N2 = integer(), Z = double(), Pval = double())
	} else {
		results = data.frame(Case1 = character(), Case2 = character(), rgname1 = character(), rgname2 = character(), Corr1 = double(), Corr2 = double(), N1 = integer(), N2 = integer(), Z = double(), Pval =double(), Qval = double())
	}
	colnames(results)[3:4] = c(rgname1, rgname2)
} else { 
	results = results[results$Pval < pcut,,drop = F]
}
write.table(pretty.numbers(results, list(Corr1..Corr2..Z = '%.3f', Pval..Qval = '%.2E')), 
	outfile, col.names = T, row.names = F, sep = "\t", quote = F)

# indata:
#    S1   S2   S3 ...
# A  1    2    1
# B  3    9    8
# ...
if (doplot && nrow(results)>0) {
	for (i in 1:nrow(results)) {
		row = results[i,]
		case1 = as.vector(unlist(row[1]))
		case2 = as.vector(unlist(row[2]))
		sam1  = as.vector(unlist(samdata[which(samdata[,2] == case1),1]))
		sam2  = as.vector(unlist(samdata[which(samdata[,2] == case2),1]))
		plotdata = as.data.frame(t(indata[unlist(row[3:4]), c(sam1,sam2),drop=F]))
		plotdata$Case = ""
		plotdata[sam1, 'Case'] = case1
		plotdata[sam2, 'Case'] = case2
		plotfile = file.path(outdir, paste(c(as.vector(unlist(row[1:4])), 'png'), collapse = '.'))
		plot.scatter(
			plotdata,
			plotfile,
			x = 1,
			y = 2,
			ggs = c(ggs, list(
				geom_smooth = list(
					aes(color = Case),
					method  = 'lm',
					se      = F
				),
				scale_color_manual = list(
					values = scales::hue_pal()(2), 
					limit  = c(case1, case2),
					labels = c(
						paste0(case1, ' (r=', round(row[5], 3),')'), 
						paste0(case2, ' (r=', round(row[6], 3),')'))
				),
				guides = list(shape = F)
			)),
			params = list(aes(shape = Case, color = Case))
		)
	}
}
