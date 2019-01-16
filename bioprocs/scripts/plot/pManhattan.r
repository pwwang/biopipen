{{rimport}}('__init__.r', 'plot.r')
options(stringsAsFactors = FALSE)

infile  = {{ i.infile | quote}}
hifile  = {{ i.hifile | quote}}
outfile = {{ o.outfile | quote}}
inopts  = {{ args.inopts | R}}
devpars = {{ args.devpars | R}}
ggs     = {{ args.ggs | R}}

# A bed file with last column the p-value and region (8 columns)
# Use region to color the snps if the snps are from the same chromosome
# chr, start, end, name, score, strand, pvalue, region
indata = read.table.inopts(infile, inopts)
pdata = data.frame(Snp = indata[,4], Chr = indata[,1], Pos = indata[,3], P = indata[,7])
if (ncol(indata) > 7)
	pdata$Region = indata[,8]

hidata = list()
if (hifile != '') {
	# possible # columns:
	# 1: snp
	# 2: snp, color
	# 3: chr, start, end
	# 4: chr, start, end, color
	hidata = read.table(hifile, header = FALSE, row.names = NULL, sep = "\t", quote = '"', check.names = FALSE)
	hincol = ncol(hidata)
	if (hincol == 1) {
		colnames(hidata) = 'snps'
	} else if (hincol == 2) {
		colnames(hidata) = c('snps', 'color')
	} else if (hincol == 3) {
		colnames(hidata) = c('chr', 'start', 'end')
	} else {
		colnames(hidata)[1:4] = c('chr', 'start', 'end', 'color')
	}
	hidata = apply(hidata, 1, as.list)
}

plot.man(pdata, plotfile = outfile, hilights = hidata, ggs = ggs, devpars = devpars)