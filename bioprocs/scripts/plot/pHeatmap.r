{{plotHeatmap}}
ggs      = {{args.ggs | Rlist}}

print ("Reading data ...")
mat      = read.table ("{{in.infile}}", sep="\t", header = {{args.header | R}}, row.names = {{args.rownames | R}}, check.names = F)
mat      = mat[order(rowSums(mat), decreasing = T), order(colSums(mat), decreasing = T)]
ssrows   = unlist(strsplit({{args.rows | quote}}, ':', fixed = T))
sscols   = unlist(strsplit({{args.cols | quote}}, ':', fixed = T))
nrows    = nrow(mat)
ncols    = ncol(mat)
dftrowN  = min(nrows, 100)
dftcolN  = min(ncols, 100)
if (length(ssrows) < 2) ssrows[2] = if (grepl("both", ssrows[1])) dftrowN/2 else dftrowN
if (length(sscols) < 2) sscols[2] = if (grepl("both", sscols[1])) dftcolN/2 else dftcolN
ssrows[2] = as.numeric(ssrows[2])
sscols[2] = as.numeric(sscols[2])

if (ssrows[1] == 'top')         rsels = 1:ssrows[2] else
if (ssrows[1] == 'bottom')      rsels = (nrows - ssrows[2]):nrows else
if (ssrows[1] == 'both')        rsels = c(1:ssrows[2], (nrows - ssrows[2]):nrows) else
if (ssrows[1] == 'random')      rsels = sample(1:nrows, ssrows[2]) else
if (ssrows[1] == 'random-both') rsels = c(sample(1:(nrows/2), ssrows[2]/2), sample((nrows - ssrows[2]):nrows, ssrows[2]/2))

if (sscols[1] == 'top')         csels = 1:sscols[2] else
if (sscols[1] == 'bottom')      csels = (ncols - sscols[2]):ncols else
if (sscols[1] == 'both')        csels = c(1:sscols[2], (ncols - sscols[2]):ncols) else
if (sscols[1] == 'random')      csels = sample(1:ncols, sscols[2]) else
if (sscols[1] == 'random-both') csels = c(sample(1:(ncols/2), sscols[2]/2), sample((ncols - sscols[2]):ncols, sscols[2]/2))

if (ssrows[1] == 'all' && sscols[1] != 'all') {
	mat = mat[, csels, drop=F]
} else if (ssrows[1] != 'all' && sscols[1] == 'all') {
	mat = mat[rsels, , drop=F]
} else if (ssrows[1] != 'all' && sscols[1] != 'all') {
	mat = mat[rsels, csels, drop=F]
}

print ("Plotting data ...")
dendro = list(dendro=TRUE, rows=NULL, cols=NULL)
{% if args.dendro | lambda x: 'dendro' in x %}
dendro$dendro = {{args.dendro['dendro'] | R}}
{% endif %}
{% if args.dendro | lambda x: 'rows' in x %}
dendro$rows = {{args.dendro['rows'] | Rvec}}
{% endif %}
{% if args.dendro | lambda x: 'cols' in x %}
dendro$cols = {{args.dendro['cols'] | Rvec}}
{% endif %}

plotHeatmap(mat, filename = {{out.outfile | quote}}, devpars = {{args.devpars | Rlist}}, ggs = {{args.ggs | Rlist}}, dendro = dendro$dendro, rows = dendro$rows, cols = dendro$cols)