from pyppl import Proc
from .utils import txt

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
pSort = Proc()
pSort.input  = "infile:file"
pSort.output = "outfile:file:{{in.infile | bn}}"
pSort.args   = {"params": "", "skip": 0}
pSort.script = """
#!/usr/bin/env bash
if [[ "{{args.skip}}" == "0" ]]; then
	sort {{args.params}} "{{in.infile}}" > "{{out.outfile}}"
else
	nhead={{args.skip}}
	ntail=$((nhead+1))
	( head -n $nhead "{{in.infile}}" && tail -n +$ntail "{{in.infile}}" | sort {{args.params}} ) > "{{out.outfile}}"
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
pFiles2Dir = Proc()
pFiles2Dir.input  = "infiles:files"
pFiles2Dir.output = "outdir:dir:{{in.infiles | [0] | fn}}.etc_{{job.index}}"
pFiles2Dir.lang   = "python"
pFiles2Dir.script = """
from glob import glob
from os import path, symlink
from sys import stderr

for fname in {{in.infiles | json}}:
	bn  = path.basename (fname)
	dst = path.join ("{{out.outdir}}", bn)
	if path.exists (dst):
		fn, _, ext = bn.rpartition('.')
		dsts = glob (path.join("{{out.outdir}}", fn + "[[]*[]]." + ext))
		print dsts
		if not dsts:
			dst2 = path.join("{{out.outdir}}", fn + "[1]." + ext)
		else:
			maxidx = max([int(path.basename(d)[len(fn)+1 : -len(ext)-2]) for d in dsts])
			dst2 = path.join("{{out.outdir}}", fn + "[" + str(maxidx+1) + "]." + ext)
		stderr.write ("Warning: rename %s to %s\\n" % (dst, dst2))
		dst = dst2
	symlink (fname, dst)
"""

"""
@name:
	pFiles2List
@description:
	Put files to a list file
@input:
	`infiles:files`: The input files
@args:
	`delimit`: The delimit. Default: r"\n"
@output:
	`outfile:file`:  The output list file
"""
pFiles2List              = Proc(desc = 'Put files to a list file.')
pFiles2List.input        = "infiles:files"
pFiles2List.output       = "outfile:file:{{in.infiles | [0] | fn}}.etc_{{job.index}}.list"
pFiles2List.args.delimit = r"\n" # r is important
pFiles2List.lang         = "python"
pFiles2List.script       = """
with open ("{{out.outfile}}", "w") as fout:
	fout.write ("{{args.delimit}}".join({{in.infiles | json}}))
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
pPat2Dir = Proc(desc = 'Link files matched by a pattern to a directory.')
pPat2Dir.input  = "pattern:var"
pPat2Dir.output = "outdir:dir:{{in.pattern | lambda x: __import__('glob').glob(x)[0] | fn }}_etc"
pPat2Dir.script = """
for fn in {{in.pattern}}; do 
	ln -s "$fn" {{out.outdir}}/
done
"""

"""
@name:
	pMergeFiles
@description:
	Merge files in the input directory
@input:
	`indir:file`: The input directory
@output:
	`outfile:file`: The output file
"""
pMergeFiles = Proc()
pMergeFiles.input  = "indir:file"
pMergeFiles.output = "outfile:file:{{indir | fn}}.merged"
pMergeFiles.script = """
> "{{out.outfile}}"
for infile in "{{indir}}/*"; do
	cat $infile >> "{{out.outfile}}"
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
pCbindList = Proc()
pCbindList.input  = "indir:file"
pCbindList.output = "outfile:file:{{indir | fn}}.mat_{{job.index}}.txt"
pCbindList.args.header = False
pCbindList.args.na     = 0
pCbindList.args.pattern= '*'
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
for (fn in Sys.glob({{args.pattern | quote}})) {
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
write.table (data, "{{out.outfile}}", col.names=T, row.names=T, sep="\\t", quote=F)
"""

pTxtFilter = Proc(desc = 'Filter a txt(tsv) file by columns and rows.')
pTxtFilter.input             = "txtfile:file"
pTxtFilter.output            = "outfile:file:{{in.txtfile | bn}}"
pTxtFilter.lang              = 'python'
pTxtFilter.args.cols         = []
pTxtFilter.args.rfilter      = 'lambda x:True'
pTxtFilter.args.header       = True
pTxtFilter.args.skip         = 0
pTxtFilter.args.delimit      = "\t"
pTxtFilter.tplenvs.txtFilter = txt.filter.py
pTxtFilter.script            = """
{{txtFilter}}

txtFilter({{in.txtfile | quote}}, {{out.outfile | quote}}, {{args.cols | json}}, {{args.rfilter}}, {{ args.header }}, {{args.skip}}, {{args.delimit | quote}})
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
pFile2Proc = Proc(desc="Convert a file to a proc so it can be used as dependent")
pFile2Proc.input  = "infile:file"
pFile2Proc.output = "outfile:file:{{in.infile | bn}}"
pFile2Proc.script = 'ln -s "{{in.infile}}" "{{out.outfile}}"'
