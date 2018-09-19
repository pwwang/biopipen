# /**
# * Calculate p-value of hypergeometric test
# * $k : overlop number
# * $m : set 1 size
# * $n : set 2 size
# * $N : background size
# *               C(m,k) * C(N-m, n-k)
# * log(P) = log[----------------------]
# *                   C(N, n)
# *               m!/(k!(m-k)!) * (N-m)!/(n-k)!(N-m-n+k)!)
# *        = log[------------------------------------------]
# *                          N!/(n!(N-n)!)
# *        =   logfact(m) + logfact(N-m) + logfact(n) + logfact(N-n)
# *          - logfact(k) - logfact(m-k) - logfact(n-k) - logfact(N-m-n+k) - logfact(N)
# */
{{rimport}}('__init__.r')
library(methods)

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
intype  = {{args.intype | R}}
N       = {{args.N | R}}
inopts  = update.list(list(rnames=T, cnames=T), {{args.inopts | R}})


ret    = matrix(nrow = 0, ncol = 5)
colnames(ret) = c('k', 'm', 'n', 'N', 'pval')

if (intype == 'raw') {
	rn     = if (inopts$rnames) 1 else NULL
	data   = read.table(infile, row.names = rn, header = inopts$cnames, sep = "\t", check.names = F)
	cnames = colnames(data)
	ncols  = ncol(data)
	if (is.null(N))
		N = nrow(data)
	if (ncols < 2)
		stop ('Need at least 2 columns.') 
	if (!inopts$cnames) {
		cnames = paste0("COL", 1:ncols)
		colnames(data) = cnames
	}
	for (i in 1:(ncols-1)) {
		for (j in (i+1):ncols) {
			row       = list(k = NA, m = NA, n = NA, N = N, pval = NA)
			row$m     = sum(data[, i], na.rm = T)
			row$n     = sum(data[, j], na.rm = T)
			row$k     = nrow(data[which(data[, i] + data[, j] == 2), ])
			row$pval  = phyper(row$k, row$m, row$N-row$m, row$n, lower.tail = F)
			row.names = paste(cnames[i], 'vs', cnames[j], sep = '.')
			row       = as.data.frame(row)

			rownames(row) = row.names

			ret = rbind(ret, row)
		}
	}
} else {
	rn     = if (inopts$rnames) 1 else NULL
	data   = read.table(infile, row.names = rn, header = inopts$cnames, sep = "\t", check.names = F)
	cnames = colnames(data)
	ncols  = ncol(data)
	if (!inopts$cnames) {
		cnames = c('k', 'm', 'n', 'N')[1:ncols]
		colnames(data) = cnames
	}
	if (!'N' %in% cnames && is.null(N)) {
		stop('N is not specified neither in infile nor args.N.')
	}
	if (!'k' %in% cnames || !'m' %in% cnames || !'n' %in% cnames) {
		stop('Column k, m or n is missing.')
	}
	rnames = rownames(data)
	if (!inopts$rnames) {
		rnames = paste0('ROW', 1:nrow(data))
		rownames(data) = rnames
	}
	for (rname in rnames) {
		if ('N' %in% cnames) N = data[rname, 'N']
		row      = list(k = NA, m = NA, n = NA, N = N, pval = NA)
		row$m    = data[rname, 'm']
		row$n    = data[rname, 'n']
		row$k    = data[rname, 'k']
		row$pval = phyper(row$k, row$m, row$N-row$m, row$n, lower.tail = F)
		row      = as.data.frame(row)

		rownames(row) = rname
		ret = rbind(ret, row)
	}
}

write.table(
	pretty.numbers(ret, formats = list(pval = '%.2E')), 
	outfile, quote = F, row.names = T, col.names = T, sep = "\t"
)