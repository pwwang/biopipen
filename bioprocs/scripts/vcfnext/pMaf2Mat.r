infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
mutypes = {{args.mutypes | R}}
binary  = {{args.binary | R}}
na      = {{args.na | R}}
samfn   = {{args.samfn}}

# read the maf file
maf  = read.table(
	if(endsWith(infile, '.gz')) gzfile(infile) else infile,
	header      = T,
	row.names   = NULL,
	check.names = F,
	sep         = "\t",
	quote       = ''
)

# select mutations with the give mutation types
if (!is.null(mutypes)) {
	maf = maf[which(maf$Variant_Classification %in% mutypes), , drop = F]
}

# get all samples
samples = unique(maf[, "Tumor_Sample_Barcode", drop = T])
# get all genes
genes   = unique(maf[, "Hugo_Symbol", drop = T])

# result matrix
ret = matrix(na, ncol = length(samples), nrow = length(genes))
rownames(ret) = genes
colnames(ret) = unique(samfn(samples))

for (gene in genes) {
	gsamples = maf[which(maf$Hugo_Symbol == gene), "Tumor_Sample_Barcode", drop = T]
	gsamples = samfn(gsamples)
	if (binary) {
		value2fill = matrix(1, nrow = length(gsamples))
		rownames(value2fill) = gsamples
	} else {
		value2fill = as.matrix(table(gsamples))
	}
	ret[gene, rownames(value2fill)] = value2fill
}
write.table(ret, outfile, row.names = T, col.names = T, sep = "\t", quote = F)
