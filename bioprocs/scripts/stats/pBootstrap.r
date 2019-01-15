{{rimport}}('__init__.r')
library(methods)
library(boot)
options(stringsAsFactors = FALSE)

{% python from os import path %}
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
outdir  = {{o.outdir | quote}}
prefix  = {{o.outdir | path.join: fn2(i.infile) | quote}}
inopts  = {{args.inopts | R}}
params  = {{args.params | R}}
nthread = {{args.nthread | R}}
# all: plot all statstics
# c(1,2): plot first 2
# FALSE or NULL: don't plot any
plots   = {{args.plot | R}}
devpars = {{args.devpars | R}}
n       = {{args.n | R}}
stats   = {{args.stats}}

# boot(data, statistic, R, sim = "ordinary", stype = c("i", "f", "w"), 
#      strata = rep(1,n), L = NULL, m = 0, weights = NULL, 
#      ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...,
#      parallel = c("no", "multicore", "snow"),
#      ncpus = getOption("boot.ncpus", 1L), cl = NULL)

# options
if (nthread > 1) {
	options(boot.parallel = 'multicore')
	options(boot.ncpus, nthread)
}

indata = read.table.inopts(infile, inopts)
statfunc = function(data, index) {
	d = data[index, ]
	r = stats(d)
	if (is.list(r)) r = unlist(r)
	as.vector(r)
}

params$data      = indata
params$statistic = statfunc
params$R         = n
b = do.call(boot, params)
rownames(b$t) = paste0('t', 1:nrow(b$t))
write.table(b$t, outfile, row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

if (!is.null(plots) && plots != FALSE) {
	if (plots == 'all') {
		plots = 1:ncol(b$t)
	}
	for (p in plots) {
		plotfile = paste0(prefix, '.stat', p, '.png')
		do.call(png, c(list(plotfile), devpars))
		plot(b, index = p)
		dev.off()
	}
}
