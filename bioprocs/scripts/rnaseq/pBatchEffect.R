
{{'plot.R', '__init__.R', 'sampleinfo.R' | rimport}}

exprfile <- {{i.expr | R}}
batchfile <- {{i.batch | R}}
tool <- {{args.tool | R}}
outfile <- {{o.outfile | R}}
outdir <- {{o.outdir | R}}
inopts <- {{args.inopts | R}}
params <- {{args.params | R}}
inlog <- {{args.inlog | R}}

expr <- read.table.inopts(exprfile, inopts)

sampleinfo <- SampleInfo2$new(batchfile)
samples <- sampleinfo$all.samples()
batch <- as.factor(sampleinfo$mat[, 'Batch'])
expr <- as.matrix(expr[, samples])
if (!inlog) {
	expr <- log2(expr + 1)
}

if (tool == 'combat') {
	library(sva)
	params$dat <- expr
	params$batch <- batch
	newexpr <- do.call(ComBat, params)
	if (!inlog) {
		newexpr <- 2^newexpr - 1
		newexpr[newexpr < 0] <- 0
	}
	write.table(round(newexpr, 3), outfile, col.names = TRUE, row.names = TRUE,
				sep="\t", quote=F)
} else { # leave it for further tools
	stop('Unsupported tool: combat.')
}
