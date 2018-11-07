library(methods)

infile   = {{i.infile | R}}
metafile = {{i.metafile | R}}
outdir   = {{o.outdir | R}}
plink    = {{args.plink | R}}
keeptxt  = {{args.keeptxt | R}}
prefix   = file.path(outdir, {{i.infile | fn2 | R}})
pedfile  = paste0(prefix, ".ped")
mapfile  = paste0(prefix, ".map")

# column names could be:
# FID, IID, PID, MID, Sex, Pheno
if (nchar(metafile) > 0) {
	metadata    = read.table(metafile, header = T, row.names = 1, sep = "\t", check.names = F)
	metasamples = rownames(metadata)
	metacols    = colnames(metadata)
} else {
	metasamples = c()
	metacols    = c()
}

indata = read.table(infile, header = T, row.names = 1, sep = "\t", check.names = F)

samples  = colnames(indata)
snps     = rownames(indata)
nsamples = ncol(indata)
nsnps    = nrow(indata)
pedncol  = nsnps + 6
pedcols  = c('FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno', snps)
comments = c('Use --compound-genotypes in later analysis')
comments = c(comments, paste(pedcols, collapse = "\t"))

mapmat = matrix(ncol = 4, nrow = nsnps)
pedmat = sapply(1:nsnps, function(i) {
	snpinfo  = unlist(strsplit(snps[i], "_"))
	mapmat[i, ] <<- c(
		if (startsWith(snpinfo[1], 'chr')) substring(snpinfo[1], 4) else snpinfo[1],
		if (startsWith(snpinfo[3], 'rs')) snpinfo[3] else paste(snpinfo[1], snpinfo[2], sep = "_"),
		0,
		snpinfo[2]
	)
	ret     = rep('00', nsamples)
	ret[which(indata[i, ] == 0)] = paste(rep(snpinfo[4], 2), collapse = '')
	ret[which(indata[i, ] == 1)] = paste(snpinfo[4:5], collapse = '')
	ret[which(indata[i, ] == 2)] = paste(rep(snpinfo[5], 2), collapse = '')
    ret
})
rm(indata)

write.table(mapmat, mapfile, sep = ' ', row.names = F, col.names = F, quote = F)
rm(mapmat)

pedmeta = matrix(NA, ncol = 6, nrow = nsamples)
rownames(pedmeta)   = samples
colnames(pedmeta)   = pedcols[1:6]
pedmeta [, 'FID']   = samples
pedmeta [, 'IID']   = samples
pedmeta [, 'PID']   = 0
pedmeta [, 'MID']   = 0
pedmeta [, 'Sex']   = 9
pedmeta [, 'Pheno'] = -9

if ('FID' %in% metacols) {
	pedmeta[metasamples, 'FID'] = metadata[metasamples, 'FID']
}
if ('PID' %in% metacols) {
	pedmeta[metasamples, 'PID'] = metadata[metasamples, 'PID']
}
if ('MID' %in% metacols) {
	pedmeta[metasamples, 'MID'] = metadata[metasamples, 'MID']
}
if ('Sex' %in% metacols) {
	pedmeta[metasamples, 'Sex'] = metadata[metasamples, 'Sex']
}
if ('Pheno' %in% metacols) {
	pedmeta[metasamples, 'Pheno'] = metadata[metasamples, 'Pheno']
}
pedmat = cbind(pedmeta, pedmat)
rm(pedmeta)

for (comment in comments) {
	cat(paste(c('#', comment, '\n')), file = pedfile, append = T)
}
write.table(pedmat, pedfile, sep = ' ', row.names = F, col.names = F, quote = F, append = T)

cmd = sprintf("%s --file %s --make-bed --out %s", plink, shQuote(prefix), shQuote(prefix))
system(cmd)

if (!keeptxt) {
	file.remove(pedfile)
	file.remove(mapfile)
}
