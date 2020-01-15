{{rimport}}('plot.r', '__init__.r')

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
fccut    = {{args.fccut | R}}
pcut     = {{args.pvalcut | R}}
usep     = {{args.usepval | R}}
hilights = {{args.hilights | R}}
devpars  = {{args.devpars | R}}
ggs      = {{args.ggs | R}}

indata = read.table.inopts(infile, list(rnames = TRUE, cnames = TRUE))
cnames = c()
for (cname in colnames(indata)) {
	cname_to_check = tolower(gsub("[^[:alnum:] ]", "", cname))
	if (cname_to_check %in% c("logfc", "log2fc", "logfoldchange", "log2foldchange")) {
		cnames = c(cnames, "logFC")
	} else if (usep && (cname_to_check %in% c("p", "pvalue", "pval"))) {
		cnames = c(cnames, "FDR")
	} else if (!usep && (cname_to_check %in% c("q", "qvalue", "qval", "fdr", "adjpval", "adjpvalue"))) {
		cnames = c(cnames, "FDR")
	} else {
		cnames = c(cnames, cname)
	}
}
colnames(indata) = cnames

plot.volplot(indata, outfile, fccut, pcut, ggs = ggs, devpars = devpars)
