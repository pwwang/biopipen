from collections import OrderedDict
from pyppl import Proc, Box
from .utils import txt, helpers, runcmd
from . import params

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
	`case`:   Case-sensitivity. Default: True
		- If True, will set $LANG as C
		- Otherwise, $LANG will be set as en_US.UTF-8
	`params`: The arguments used by `sort`
"""
pSort                        = Proc(desc = 'Sort file.')
pSort.input                  = "infile:file"
pSort.output                 = "outfile:file:{{in.infile | bn}}"
pSort.args.params            = OrderedDict()
pSort.args.skip              = 0
pSort.args.case              = True
pSort.tplenvs.runcmd         = runcmd.py
pSort.tplenvs.params2CmdArgs = helpers.params2CmdArgs.py
pSort.lang                   = params.python.value
pSort.script                 = "file:scripts/common/pSort.py"

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
pFiles2Dir = Proc(desc = 'Put files to a directory using symbolic links.')
pFiles2Dir.input  = "infiles:files"
pFiles2Dir.output = "outdir:dir:{{in.infiles | lambda x: sorted(x) | [0] | fn}}.etc"
pFiles2Dir.lang   = params.python.value
pFiles2Dir.script = "file:scripts/common/pFiles2Dir.py"

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
pMergeFiles              = Proc(desc = 'Merge files.')
pMergeFiles.input        = "infiles:files"
pMergeFiles.output       = "outfile:file:{{in.infiles[0] | fn}}_etc-merged{{in.infiles[0] | ext}}"
pMergeFiles.args.skip    = []
pMergeFiles.args.comment = '#'
pMergeFiles.args.header  = False
pMergeFiles.lang         = params.python.value
pMergeFiles.script       = "file:scripts/common/pMergeFiles.py"

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
pTxtFilter.lang              = params.python.value
pTxtFilter.args.cols         = []
pTxtFilter.args.rfilter      = None
pTxtFilter.args.header       = True
pTxtFilter.args.skip         = 0
pTxtFilter.args.delimit      = "\t"
pTxtFilter.tplenvs.txtFilter = txt.filter.py
pTxtFilter.script            = "file:scripts/common/pTxtFilter.py"


pTxtTransform                      = Proc(desc = 'Transform a txt(tsv) file.')
pTxtTransform.input                = "txtfile:file"
pTxtTransform.output               = "outfile:file:{{in.txtfile | bn}}"
pTxtTransform.lang                 = params.python.value
pTxtTransform.args.transformer     = None
pTxtTransform.args.header          = True
pTxtTransform.args.delimit         = "\t"
pTxtTransform.args.comment         = "#"
pTxtTransform.tplenvs.txtTransform = txt.transform.py
pTxtTransform.script               = "file:scripts/common/pTxtTransform.py"

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

"""
@name:
	pStr2File
@description:
	Save string to a file.
@input:
	`in:var`: The input string.
@output:
	`outfile:file`: The output file.
"""
pStr2File                   = Proc(desc = "Save string to a file.")
pStr2File.input             = "in:var"
pStr2File.output            = "outfile:file:{{in.in | encode}}.txt"
pStr2File.args.breakOn      = ','
pStr2File.args.trimLine     = True
pStr2File.tplenvs.encode    = lambda x: __import__('re').sub(r'[^\w_]', '', x)[:16]
pStr2File.lang              = params.python.value
pStr2File.script            = "file:scripts/common/pStr2File.py"

"""
@name:
	pSimRead
@description:
	Read files 
@input:
	`infiles:files`: The input files
@output:
	`outfile:file`: The output file
@args:
	`skip`: argument skip for each file
	`delimit`: argument delimit for each file
	`gzip`: argument gzip for each file
	`match`: The match function. 
	`do`: The do function. Global vaiable `fout` is available to write results to output file.
@requires:
	[`python-simread`](https://github.com/pwwang/simread)
"""
pSimRead              = Proc(desc = 'Read files simultaneously.')
pSimRead.input        = 'infiles:files'
pSimRead.output       = 'outfile:file:{{in.infiles[0] | fn}}.etc.simread.txt'
pSimRead.args.skip    = []
pSimRead.args.delimit = []
pSimRead.args.gzip    = []
pSimRead.args.match   = None
pSimRead.args.do      = None
pSimRead.lang         = params.python.value
pSimRead.script       = "file:scripts/common/pSimRead.py"

pAddHeader        = Proc(desc = 'Add the header of 1st file to 2nd file.')
pAddHeader.input  = "infile1:file, infile2:file"
pAddHeader.output = "outfile:file:{{in.infile2 | bn}}"
pAddHeader.args.n = 1
pAddHeader.script = "file:scripts/common/pAddHeader.bash"
