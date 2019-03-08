{{rimport}}('__init__.r')

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
outdir  = {{o.outdir | quote}}
inopts  = {{args.inopts | R}}
plotman = {{args.plot | R}}
ggs     = {{args.ggs | R}}
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
hifile  = {{args.hifile | R}}
hilabel = {{args.hilabel | R}}

library(coloc)
options(stringsAsFactors = FALSE)

indata = read.table.inopts(infile, inopts)
snps   = indata[,4]
nsnps  = nrow(indata)
ret = coloc.abf(
	dataset1 = list(pvalues = indata[, 8], snp = snps, N = nsnps, type = "quant"),
	dataset2 = list(pvalues = indata[, 9], snp = snps, N = nsnps, type = "quant"),
	MAF      = indata[, 7]
)
write.table(as.data.frame(ret$summary), outfile, col.names = FALSE, sep = "\t", quote = FALSE)
detailfile = file.path(outdir, 'coloc.details.txt')
write.table(ret$results, detailfile, sep = "\t", quote = FALSE, row.names = FALSE)

if (!plotman) 
	q(save = FALSE)

{{rimport}}('plot.r')

library(ggrepel)
# assign regions
nbins   = 10
minpos  = min(indata[,3])
maxpos  = max(indata[,3])
binsize = (maxpos - minpos) / nbins
bits    = round(log10(binsize))
nbits   = binsize/10 ** bits
binsize = round(nbits * 4) / 4 * 10 ** bits
nbins   = ceiling((maxpos - minpos)/binsize)
start   = as.integer(minpos / binsize + 1) * binsize 
glabels = seq(start, maxpos, binsize)
gbreaks = glabels - minpos
glabels = sapply(glabels, function(x) {
	if (x < 1000) {
		return(x)
	} else if (x < 1000000) {
		return(paste0(x/1000, 'K'))
	} else {
		return(paste0(x/1000000, 'M'))
	}
})

if (is.true(colnames(indata))) {
	phenotypes = colnames(indata)[8:9]
} else {
	phenotypes = c('Phenotype1', 'Phenotype2')
}
data1 = data.frame(X = indata[,3] - min(indata[,3]) + 1, Y = -log10(indata[, 8]), Group = phenotypes[1])
data2 = data.frame(X = data1$X, Y = -log10(indata[, 9]), Group = phenotypes[2])
data  = rbind(data1, data2)
rm(data1)
rm(data2)

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

ggs_hilight = list()
hinames = names(hidata)
if (!is.null(hinames)) {
	hidata = list(hidata)
}
hdata = NULL
for (hilight in hidata) {
	if ('snps' %in% names(hilight)) {
		hdata_tmp = data[Snp %in% hilight$snps]
		hdata = ifelse(is.null(hdata), hdata_tmp, rbind(hdata, hdata_tmp))
		rm(hdata_tmp)
	} else if ('chr' %in% names(hilight)) {
		hdata_tmp = data[Chr == hilight$chr & Pos >= as.numeric(hilight$start) & Pos <= as.numeric(hilight$end)]
		hdata = ifelse(is.null(hdata), hdata_tmp, rbind(hdata, hdata_tmp))
		rm(hdata_tmp)
	} else {
		stop('Either "snps" or "chr" needed to locate the snps to be highlighted.')
	}
}
if (!is.null(hdata)) {
	hicolor = ifelse(is.null(hilight$color), "red", hilight$color)
	ggs_hilabel = list()
	if (hilabel) {
		ggs_hilabel = list(geom_label_repel = list(aes_string(x = 'X', y = 'Y', label = 'Snp'), color = "white", fill = hicolor, data = hdata, inherit.aes = FALSE, alpha = .8))
	}
	ggs_hilight = c(
		ggs_hilight, 
		ggs_hilabel,
		list(geom_point = list(aes_string(x = 'X', y = 'Y'), color = hicolor, data = hdata, inherit.aes = FALSE))
	)
}
ggs = c(list(
	geom_vline         = list(xintercept = gbreaks, color = '#DEDEDE', size = .3),
	expand_limits      = list(y = c(0, max(data[, 'Y']) * 1.05)),
	geom_point         = list(aes_string(color = 'Group'), alpha = .8),
	scale_y_continuous = list(expand = c(0, 0)),
	scale_x_continuous = list(expand = expand_scale(mult = c(0.01, 0.01)), label = glabels, breaks= gbreaks),
	theme_bw           = list(),
	theme              = list(
		panel.grid.minor   = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5)
	),
	labs               = list(x = "", y = "-log10(p-value)")
	), ggs_hilight, ggs)
plotfile = file.path(outdir, 'coloc.man.png')
plot.scatter (data, plotfile, x = 'X', y = 'Y', params = list(aes_string(color = "Group")), ggs = ggs, devpars = devpars)



