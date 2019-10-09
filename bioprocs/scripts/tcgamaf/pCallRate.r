{{runcmd}}
{{cbindfill}}
{{plotHist}}

outsample = file.path("{{o.outdir}}", "{{i.indir | fn}}.sampleCallRate.txt")
figsample = file.path("{{o.outdir}}", "{{i.indir | fn}}.sampleCallRate.png")
outsnp    = file.path("{{o.outdir}}", "{{i.indir | fn}}.snpCallRate.txt")
figsnp    = file.path("{{o.outdir}}", "{{i.indir | fn}}.snpCallRate.png")

setwd("{{i.indir}}")
files = list.files(pattern = "*.vcf")

data  = NULL

for (file in files) {
	sample  = tools::file_path_sans_ext(tools::file_path_sans_ext(file))
	tmp     = runcmd (paste('awk', '\'BEGIN{OFS="\\t"} $0 !~ /^#/ {print $1":"$2,1}\'', file), intern=T)
	con     =  textConnection(tmp)
	tmp     = read.table (con, header=F, row.names=1, check.names=F)
	close(con)
	colnames(tmp) = c(sample)
	data    = cbindfill(data, tmp)
}

samplecr = colSums(data) / nrow(data)
snpcr    = rowSums(data) / ncol(data)

# sample/snp call rate
write.table (samplecr, outsample, quote=F, sep="\t", col.names=F)
plotHist (samplecr, figsample, devpars={{args.devpars | Rlist}}, ggs={{args.histplotggs | lambda x: x + ['r:ggtitle("Sample call rate")', 'r:ylab("# samples")', 'r:xlab("Call rate")'] | Rlist}})
write.table (snpcr, outsnp, quote=F, sep="\t", col.names=F)
plotHist (snpcr, figsnp, devpars={{args.devpars | Rlist}}, ggs={{args.histplotggs | lambda x: x + ['r:ggtitle("SNP call rate")', 'r:ylab("# SNPs")', 'r:xlab("Call rate")'] | Rlist}})