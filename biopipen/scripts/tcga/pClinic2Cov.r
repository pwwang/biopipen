{{rimport}}('__init__.r')

{% python from bioprocs.utils import alwaysList %}
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
covs    = {{args.covs | alwaysList | R}}
asnum   = {{args.asnum | R}}
suffix  = {{args.suffix | alwaysList | R}}

if (length(suffix) == 0) suffix = ''

cldata = read.table.inopts(infile, list(
	sep = "\t",
	skip = 1,
	cnames = TRUE,
	rnames = FALSE
))

# remove CDE_ID: row
covdata = cldata[-1, covs, drop = FALSE]
covdata[covdata == '']                = '[Not Applicable]'
covdata[covdata == '[Not Available]'] = '[Not Applicable]'
covdata[covdata == '[Not Evaluated]'] = '[Not Applicable]'
covdata[covdata == '[Unknown]']       = '[Not Applicable]'
if (asnum) {
	for (cov in covs) {
		covdata[, cov] = as.numeric(covdata[, cov])
	}
}

ret = NULL
for (suf in suffix) {
	rownames(covdata) = paste0(cldata[-1, 'bcr_patient_barcode'], suf)
	if (is.null(ret)) {
		ret = covdata
	} else {
		ret = rbind(ret, covdata)
	}
}

write.table(ret, outfile, sep = "\t", quote = FALSE)
