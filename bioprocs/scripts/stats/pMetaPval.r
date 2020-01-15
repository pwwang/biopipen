library("methods")
library("metap")
{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
intype  = {{args.intype | R}}
inopts  = {{args.inopts | R}}
method  = {{args.method | R}}
na      = {{args.na | R}}
if (method == 'fisher') {
	method = 'sumlog'
}

indata = read.table.inopts(infile, inopts)
if (intype == 'matrix') {
	entries = rownames(indata)
	if (is.null(entries)) {
		entries = 1:nrow(indata)
	}
} else {
	entries = levels(factor(indata[,1]))
}

ret = data.frame(Entry = character(), N = numeric(), MetaPval = numeric())
for (entry in entries) {
	if (intype == 'matrix') {
		pvals = indata[entry,]
	} else {
		pvals = indata[which(indata[,1] == entry), 2]
	}
	pvals = na.omit(pvals)
	lenps = length(pvals)
	if (lenps < 1) {
		ret = rbind(ret, list(Entry = entry, N = 0, MetaPval = NA))
	} else if (lenps == 1) {
		ret = rbind(ret, list(Entry = entry, N = 1, MetaPval = pvals))
	} else {
		ret = rbind(ret, list(Entry = entry, N = length(pvals), MetaPval = do.call(method, list(pvals))$p))
	}
}

write.table(pretty.numbers(ret, list(MetaPval = '%.2E')), outfile, row.names = F, col.names = T, sep = "\t", quote = F)
