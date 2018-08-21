infile  = {{in.infile | quote}}
outfile = {{out.outfile | quote}}
style   = {{args.style | R}}
inopts  = {{args.inopts | R}}
acode   = {{args.acode | R}}
bcode   = {{args.bcode | R}}

bcode = c("# Before table inserted:", bcode)
acode = c("# After table inserted:", acode)

cat (bcode, file = outfile, sep = "\n")

rtopts         = inopts
rtopts$delimit = NULL
rtopts$rnames  = NULL
rtopts$cnames  = NULL

rtopts$sep = inopts$delimit
rtopts$header    = inopts$cnames
if (inopts$rnames) {
	rtopts$row.names = 1
} else {
	rtopts = c(rtopts, list(row.names = NULL))
}
rtopts = c(infile, rtopts)
mat    = do.call(read.table, rtopts)
nrows  = nrow(mat)
ncols  = ncol(mat)

doctablecols = ncols
if (inopts$rnames) doctablecols = ncols + 1
doctablerows = nrows
if (inopts$cnames) doctablerows = nrows + 1

codes  = c("# Inserting table")
codes  = c(codes, sprintf("table = doc.add_table(rows = %d, cols = %d)", doctablerows, doctablecols))
if (!is.null(style)) {
	codes  = c(codes, sprintf("table.style = '%s'", style))
}
if (inopts$cnames) {
	cnames = colnames(mat)
	codes = c(codes, "header = table.rows[0].cells")
	if (inopts$rnames) {
		codes = c(codes, "header[0].text = ''")
		for (i in 1:ncols) {
			codes = c(codes, sprintf("header[%s].text = '''%s'''", i, cnames[i]))
		}
	} else {
		for (i in 1:ncols) {
			codes = c(codes, sprintf("header[%s].text = '''%s'''", i-1, cnames[i]))
		}
	}
}
rnames = rownames(mat)
for (i in 1:nrows) {
	rindex = if (inopts$cnames) i else i-1
	rname  = rnames[i]
	rdata  = mat[i,,drop=T]
	codes = c(codes, sprintf("row = table.rows[%d].cells", rindex))
	if (inopts$rnames) {
		codes = c(codes, sprintf("row[0].text = '''%s'''", rname))
		for (j in 1:ncols) {
			codes = c(codes, sprintf("row[%d].text = '''%s'''", j, unlist(rdata[j])))
		}
	} else {
		for (j in 1:ncols) {
			codes = c(codes, sprintf("row[%d].text = '''%s'''", j-1, unlist(rdata[j])))
		}
	}
}

cat (codes, file = outfile, sep = "\n", append = T)
cat (acode, file = outfile, sep = "\n", append = T)
