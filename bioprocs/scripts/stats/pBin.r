
{{rimport}}('__init__.r')
{% assign default_inopts = {'delimit': '\t', 'rnames': False, 'cnames': True, 'check.names': False} %}
{% python default_inopts.update(args.inopts) %}
{% python default_inopts['sep'] = default_inopts['delimit'] %}
{% python default_inopts['row.names'] = default_inopts['rnames'] or None %}
{% python default_inopts['header'] = default_inopts['cnames'] %}
{% python del default_inopts['delimit'] %}
{% python del default_inopts['rnames'] %}
{% python del default_inopts['cnames'] %}
{% assign default_binopts = {'nbin': None, 'step': None, 'nan': 'skip', 'out': 'lower'} %}
{% python default_binopts.update(args.binopts) %}

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{default_inopts | R}}
binopts = {{default_binopts | R}}
cols    = {{args.cols | R: False}}

inopts = c(infile, inopts)
data   = do.call(read.table, inopts)
cnames = colnames(data)

newcols = list()
if (is.null(cols)) {
	for (col in cnames) {
		newcols[[col]] = binopts
	}
} else if (is.list(cols)) {
	for (col in names(cols)) {
		cname = as.integer(col)
		if (is.na(cname)) 
			cname = col
		else
			cname = cnames[cname]
		newcols[[cname]] = update.list(binopts, cols[[col]])
	}
} else {
	for (col in cols) {
		cname = as.integer(col)
		if (is.na(cname)) 
			cname = col
		else
			cname = cnames[cname]
		newcols[[cname]] = binopts
	}
}

ret = NULL
for (i in 1:ncol(data)) {
	cname = cnames[i]
	dcol  = data[, i, drop = F]
	if (cname %in% names(newcols)) {
		attach (newcols[[cname]])
		if (nan == 'as0') {
			dcol[is.na(dcol)] = 0
			dcol2 = dcol
		} else {
			dcol2 = dcol[!is.na(dcol)]
		}
		
		cmax = max(dcol2)
		cmin = min(dcol2)
		probs = NULL
		if (!is.null(nbin)) {
			breaks = quantile(dcol2, probs = seq(0, 1, 1/nbin))
		} else if (!is.null(step)) {
			breaks = seq(cmin, cmax + 1, step)
		} else {
			stop('One of `nbin` and `step` is required for binning.')
		} 
		cuts = cut(dcol2, breaks = breaks, include.lowest = T)
		ucuts = levels(factor(cuts))
		for (j in 1:length(ucuts)) {
			idx = which(cuts == ucuts[j])
			if (out == 'step') {
				dcol2[idx] = unlist(strsplit(substring(ucuts[j], 2), ',', fixed=T))[1]
			} else if (out == 'lower' || out == 'min') {
				dcol2[idx] = min(dcol2[idx])
			} else if (out == 'upper' || out == 'max') {
				dcol2[idx] = max(dcol2[idx])
			} else if (out == 'mean') {
				dcol2[idx] = mean(dcol2[idx])
			} else if (out == 'median') {
				dcol2[idx] = median(dcol2[idx])
			} else {
				dcol2[idx] = j
			}
		}
		dcol[!is.na(dcol)] = dcol2
	}
	if (is.null(ret)) {
		ret = dcol
	} else {
		ret = cbind(ret, dcol)
	}
}

write.table(ret, outfile, sep = inopts$sep, col.names = inopts$header, row.names = if (is.null(inopts$row.names)) F else T, quote = F)




