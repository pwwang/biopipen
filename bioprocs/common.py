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
pSort.output = "outfile:file:{{infile | bn}}.sorted"
pSort.args   = {"params": "", "skip": 0}
pSort.script = """
#!/usr/bin/env bash
if [[ "{{args.skip}}" == "0" ]]; then
	sort {{args.params}} "{{infile}}" > "{{outfile}}"
else
	nhead={{args.skip}}
	ntail=$((nhead+1))
	( head -n $nhead "{{infile}}" && tail -n +$ntail "{{infile}}" | sort {{args.params}} ) > "{{outfile}}"
fi
"""

"""
@name:
	pFiles2Dir
@description:
	A helper process to convert a list of files into a directory, so that some processes can take it as input
@input:
	`infiles:files`: The input files
@output:
	`outdir:dir`:    The output directory
"""
pFiles2Dir = proc()
pFiles2Dir.input  = "infiles:files"
pFiles2Dir.output = "outdir:dir:{{infiles | [0] | fn}}.etc.{{#}}"
pFiles2Dir.lang   = "python"
pFiles2Dir.script = """
from glob import glob
from os import path, symlink
from sys import stderr

for fn in [{{infiles | acquote}}]:
	bn  = path.basename (fn)
	dst = path.join ("{{outdir}}", bn)
	if path.exists (dst):
		dsts = glob (path.join("{{outdir}}", bn + "[*]"))
		if not dsts:
			dst2 = path.join("{{outdir}}", bn + "[1]")
		else:
			maxIdx = max ([int(path.basename(d)[:-1].rpartition('[')[-1]) for d in dsts])
			dst2 = path.join("{{outdir}}", bn + "[%s]" % str(maxIdx+1))
		stderr.write ("Warning: rename %s to %s" % (dst, dst2))
		dst = dst2
	symlink (fn, dst)
"""

"""
@name:
	pPat2Dir
@description:
	A helper process to convert a list of files by a pattern (wildcards) into a directory, so that some processes can take it as input
@input:
	`pattern:var`: The pattern
@output:
	`outdir:dir`:    The output directory
"""
pPat2Dir = proc()
pPat2Dir.input  = "pattern:var"
pPat2Dir.output = "outdir:dir:{{pattern | __import__('glob').glob(_) | [0] | fn }}.etc.{{#}}"
pPat2Dir.script = """
for fn in {{pattern}}; do
	ln -s "$fn" "{{outdir}}/"
done
"""

"""
@name:
	pMergeFile
@description:
	Merge files in the input directory
@input:
	`indir:file`: The input directory
@output:
	`outfile:file`: The output file
"""
pMergeFile = proc()
pMergeFile.input  = "indir:file"
pMergeFile.output = "outfile:file:{{indir | fn}}.merged"
pMergeFile.script = """
> "{{outfile}}"
for infile in "{{indir}}/*"; do
	cat $infile >> "{{outfile}}"
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
pCbindList.output = "outfile:file:{{indir | fn}}.{{#}}.mat.txt"
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
header = {{args.header | Rbool}}
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
if (!is.na({{args.na}})) {
	data[is.na(data)] = {{args.na}}
}
write.table (data, "{{outfile}}", col.names=T, row.names=T, sep="\\t", quote=F)
"""

"""
@name:
	pFile2Proc
@description:
	Convert a file to a proc so it can be used as dependent
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
"""
pFile2Proc = proc (desc="Convert a file to a proc so it can be used as dependent")
pFile2Proc.input  = "infile:file"
pFile2Proc.output = "outfile:file:{{infile | bn}}"
pFile2Proc.script = 'ln -s "{{infile}}" "{{outfile}}"'
