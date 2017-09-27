from pyppl import Proc


library(simpleaffy)
setwd ("/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pCeldir2Matrix.notag.2XgyTuAe/0/input/17158004339_1O-08miRNA_etc")
files    = list.celfiles()

covdescfile = file.path("/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pCeldir2Matrix.notag.2XgyTuAe/0/input/17158004339_1O-08miRNA_etc", "covdesc.txt")
# generate a covdesc file
if ("" == '') {
  covdesc  = matrix(1, ncol = 1, nrow = length(files))
  colnames(covdesc) = c("treatment")
  rownames(covdesc) = files
  write.table(covdesc, covdescfile, quote = F, sep = "\t")
} else {
  file.link("", covdescfile)
}

affydata = read.affy(covdesc = "covdesc.txt")
if ("/data2/junwenwang/panwen/deFilipisMiRNA/raw/annotations/miRNA-4_0-st-v1.cdf" != '') {
  library(makecdfenv)
  cdffile = "/data2/junwenwang/panwen/deFilipisMiRNA/raw/annotations/miRNA-4_0-st-v1.cdf"
  cdfname = cleancdfname(whatcdf(files[1]))
  assign(cdfname, make.cdf.env(basename(cdffile), cdf.path = dirname(cdffile)))
  affydata@cdfName = cdfname
} 
exprs    = call.exprs(affydata, algorithm = "rma")
exprsout = exprs@assayData$exprs

# annotate

annos    = read.table("/data2/junwenwang/panwen/deFilipisMiRNA/raw/annotations/miRNA-4_0-st-v1.annotations.20160922.csv", header=T, sep=",",colClasses = c("NULL", NA, "NULL", NA, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL"), row.names = 2, check.names = F )
rnames   = intersect(rownames(exprsout), rownames(annos))
exprsout = exprsout[rnames,,drop=F]
rownames(exprsout) = make.unique(as.vector(annos[rnames,,drop=T]))


cnames   = colnames(exprsout)
cnames   = gsub("\\.CEL$", "", cnames)
cnames   = gsub("\\.cel$", "", cnames)
colnames(exprsout) = cnames
write.table(exprsout, "/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pCeldir2Matrix.notag.2XgyTuAe/0/output/17158004339_1O-08miRNA_etc/17158004339_1O-08miRNA_etc.rma.expmat", quote=F, sep="\t")

# do qc


####


# require sva to remove batch effect
# get the exp data
ematrix    = read.table("/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pDeg.notag.usDyRdCi/0/input/17158004339_1O-08miRNA_etc.rma.expmat",  header=T, row.names = 1, check.names=F, sep="\t")
filters    = c(5,3)

# get group names
groupmat   = read.table("/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pDeg.notag.usDyRdCi/0/input/17158004339_1O-08miRNA_etc.txt", header = TRUE, check.names = F, sep = "\t", row.names = 1)
groups     = unique(groupmat[, 1, drop=T])

datas      = list()
for (group in groups) {
	samples = rownames(groupmat[groupmat == group,,drop=F])
	samples = gsub('\\.CEL$', '', samples, ignore.case = T)
	datas   = c(datas, list(ematrix[, samples]))
}
rm(ematrix)



for (i in 1:length(groups)) {
	for (j in 1:length(groups)) {
		if (j <= i) next
		# paths
		outdir  = file.path("/data2/junwenwang/panwen/deFilipisMiRNA/workdir/PyPPL.pDeg.notag.usDyRdCi/0/output/17158004339_1O-08miRNA_etc.rma-17158004339_1O-08miRNA_etc", paste(groups[i], groups[j], sep='-'))
		allfile = file.path(outdir, '17158004339_1O-08miRNA_etc.rma-17158004339_1O-08miRNA_etc.all')
		outfile = file.path(outdir, '17158004339_1O-08miRNA_etc.rma-17158004339_1O-08miRNA_etc.degs')
		vcplot  = file.path(outdir, '17158004339_1O-08miRNA_etc.rma-17158004339_1O-08miRNA_etc.volcano.png')
		hmplot  = file.path(outdir, '17158004339_1O-08miRNA_etc.rma-17158004339_1O-08miRNA_etc.heatmap.png')
		dir.create(outdir, showWarnings = F, recursive = T)

		
		library(limma)
		# design
		design = model.matrix(~factor(groupmat[, 1, drop=T]))
		colnames(design) = c(groups)

		# fit
		data   = cbind(datas[[i]], datas[[j]])
		data   = data[rowSums(data>filters[1]) >= filters[2],,drop=F]
		fit    = lmFit(data, design)
		fit    = eBayes(fit)

		# degs
		alls   = topTable(fit, coef=1, n=nrow(data), adjust="BH")
		degs   = alls[alls$adj.P.Value<0.05,,drop=F]
		write.table(alls, allfile, quote=F, sep="	")
		write.table(degs, outfile, quote=F, sep="	")
		rm(alls)
		
		# volcano plot
		
		png (vcplot, res=300, width=2000, height=2000)
		volcanoplot(fit, highlight = nrow(degs))
		dev.off()
		

		# heatmap
		
		## indent: keep ##
		
require('pheatmap')
if (!exists('plotHeatmap')) {
	plotHeatmap = function (mat, params = list()) {
		do.call(pheatmap, c(list(mat=mat), params))
	}
}

		tmatrix    = datas[[i]][rownames(degs),,drop=F]
		nmatrix    = datas[[j]][rownames(degs),,drop=F]
		tmatrix    = tmatrix + 1
		nmatrix    = rowSums(nmatrix)/ncol(nmatrix)
		nmatrix    = nmatrix + 1
		log2fc     = log2(apply(tmatrix, 2, function(col) col/nmatrix))
		plotHeatmap(log2fc, list(filename=hmplot, show_rownames=FALSE))
		

		
	}
}⏎