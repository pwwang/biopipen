{{'__init__.r' | rimport}}

infile = {{i.infile | R}}
outdir = {{o.outdir | R}}
inopts = {{args.inopts | R}}
by     = {{args.by}}
fntemp = {{o.outdir, stem2(i.infile) | @join: "/"
									 | @append: ".<placeholder>"
									 | @append: ext(i.infile) | quote}}

mat    = read.table.inopts(infile, inopts)
cnames = colnames(mat)

if (is.numeric(by)) {
	if (is.null(cnames)) {
		cnames = paste0("COL", 1:ncol(mat))
		colnames(mat) = cnames
	}
	for ( chunk in split(1:length(cnames), ceiling(seq_along(1:length(cnames))/by)) ) {
		cname = cnames[chunk]
		m     = mat[, cname, drop = FALSE]
		cname = paste(gsub("[[:punct:]]", "_", cname), collapse='-')
		cname = gsub("_+", "_", cname)
		fn    = gsub("<placeholder>", cname, fntemp)
		write.table(m, fn, quote = FALSE, sep = "\t", row.names = inopts$rnames, col.names = inopts$cnames)
	}
} else if (is.function(by)) {
	groups = list()
	if (is.null(cnames)) {
		cnames = 1:ncol(mat)
	}
	for (cname in cnames) {
		if (is.null(groups[[by(cname)]])) {
			groups[[by(cname)]] = cname
		} else {
			groups[[by(cname)]] = c(groups[[by(cname)]], cname)
		}
	}
	for (group in names(groups)) {
		cname = groups[[group]]
		m     = mat[, cname, drop = FALSE]
		cname = paste(gsub("[[:punct:]]", "_", cname), collapse='-')
		cname = gsub("_+", "_", cname)
		fn    = gsub("<placeholder>", cname, fntemp)
		write.table(m, fn, quote = FALSE, sep = "\t", row.names = inopts$rnames, col.names = inopts$cnames)
	}
} else {
	stop('Unsupported `args.by`.')
}
