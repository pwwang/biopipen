from pyppl import proc

"""
@name:
	pSort
@description:
	Sort file using linux command `sort`
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`skip`:   To skip first N lines. Default: 0
	`params`: The arguments used by `sort`
"""
pSort = proc ()
pSort.input  = "infile:file"
pSort.output = "outfile:file:{{infile.bn}}.sorted"
pSort.args   = {"params": "", "skip": 0}
pSort.script = """
#!/usr/bin/env bash
if [[ "{{proc.args.skip}}" == "0" ]]; then
	sort {{proc.args.params}} "{{infile}}" > "{{outfile}}"
else
	nhead={{proc.args.skip}}
	ntail=$((nhead+1))
	( head -n $nhead "{{infile}}" && tail -n +$ntail "{{infile}}" | sort {{proc.args.params}} ) > "{{outfile}}"
fi
"""

"""
@name:
	pFiles2Dir
@description:
	A helper process to convert a list of file into a directory, so that some processes can take it as input
@input:
	`infiles:files`: The input files
@output:
	`outdir:dir`:    The output directory
"""
pFiles2Dir = proc()
pFiles2Dir.input  = "infiles:files"
pFiles2Dir.output = "outdir:dir:{{infiles.fn | [0]}}_etc.{{#}}"
pFiles2Dir.script = """
for fn in "{{infiles | ('" "').join(_)}}"; do
	ln -s "$fn" "{{outdir}}/"
done
"""

"""
@name:
	pCbindList
@description:
	Column bind lists, fill miss rows with specific value
@input:
	`indir:file`: The directory containing the list files
	- header can be omited, but row names are required
@output:
	`outfile:file`: The output matrix
@args:
	`header`: Whether list has header. Default: False (will use file name as header)
	`na`:     The missing values. Default: 0
	- If it's a string, remember the quote (i.e.: '"missing"')
"""
pCbindList = proc ()
pCbindList.input  = "indir:file"
pCbindList.output = "outfile:file:{{indir.fn}}.{{#}}.mat.txt"
pCbindList.args   = {"header": False, "na": "0"}
pCbindList.lang   = "Rscript"
pCbindList.script = """
cbind.fill = function (x1, x2) {
    y = merge(x1, x2, by='row.names', all=T, sort=F)
    rownames(y) = y[, "Row.names"]
    y = y[, -1, drop=F]
    cnames      = c(colnames(x1), colnames(x2))
    if (!is.null(cnames)) {
        colnames(y) = cnames
    }
    return (y)
}

setwd ("{{indir}}")
data   = NULL
header = {{proc.args.header | str(_).upper()}}
for (fn in list.files()) {
	x = read.table (fn, header = header, row.names=1, check.names=F, sep="\\t")
	if (!header) {
		colnames(x) = c(unlist(strsplit(fn, ".", fixed=T))[1])
	}
	if (is.null(data)) {
		data = x
	} else {
		data = cbind.fill (data, x)
	}
}
if (!is.na({{proc.args.na}})) {
	data[is.na(data)] = {{proc.args.na}}
}
write.table (data, "{{outfile}}", col.names=T, row.names=T, sep="\\t", quote=F)
"""