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
	`mem`    : The buffer size. Default: 4G
	`tmpdir` : The tmpdir.
	`unique` : Just keep the unique lines. Default: False
	`delimit`: The delimit to separate the fields. Default: '\t'
	`noeline`: Remove empty lines. Default: False
	`params` : The arguments used by `sort`
"""
pSort                     = Proc(desc = 'Sort file.')
pSort.input               = "infile:file"
pSort.output              = "outfile:file:{{in.infile | bn}}"
pSort.args.params         = OrderedDict()
pSort.args.skip           = 0
pSort.args.case           = True
pSort.args.mem            = params.mem4G.value
pSort.args.tmpdir         = params.tmpdir.value
pSort.args.unique         = False
pSort.args.delimit        = "\t"
pSort.args.noeline        = False
pSort.envs.runcmd         = runcmd.py
pSort.envs.params2CmdArgs = helpers.params2CmdArgs.py
pSort.lang                = params.python.value
pSort.script              = "file:scripts/common/pSort.py"

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
pFiles2Dir                = Proc(desc = 'Put files to a directory using symbolic links.')
pFiles2Dir.input          = "infiles:files"
pFiles2Dir.output         = "outdir:dir:{{in.infiles | lambda x: sorted(x) | [0] | fn}}.dir"
pFiles2Dir.lang           = params.python.value
pFiles2Dir.script         = "file:scripts/common/pFiles2Dir.py"

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
pFile2Proc                = Proc(desc="Convert a file to a proc so it can be used as dependent")
pFile2Proc.input          = "infile:file"
pFile2Proc.output         = "outfile:file:{{in.infile | bn}}"
pFile2Proc.script         = 'ln -s "{{in.infile}}" "{{out.outfile}}"'

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
pStr2File                 = Proc(desc = "Save string to a file.")
pStr2File.input           = "in:var"
pStr2File.output          = "outfile:file:{{in.in | encode}}.txt"
pStr2File.args.breakOn    = ','
pStr2File.args.trimLine   = True
pStr2File.envs.encode     = lambda x: __import__('re').sub(r'[^\w_]', '', x)[:16]
pStr2File.lang            = params.python.value
pStr2File.script          = "file:scripts/common/pStr2File.py"

"""
@name:
	pPrepend
@description:
	Prepend a string to a file
@input:
	`in:var`: The input string.
	`infile:file`: The input file.
@output:
	`outfile:file`: The output file.
"""
pPrepend                  = Proc(desc = "Save string to a file.")
pPrepend.input            = "in:var, infile:file"
pPrepend.output           = "outfile:file:{{in.infile | fn}}.prepend.txt"
pPrepend.script           = "file:scripts/common/pPrepend.bash"

"""
@name:
	pAddHeader
@description:
	Add the header of 1st file to 2nd file.
@input:
	`infile1:file`: The first file containing the header.
	`infile2:file`: The second file with the body.
@output:
	`outfile:file`: The output file with the header from 1st input file, body from 2nd file.
@args:
	`n`: The number of header lines.
"""
pAddHeader                = Proc(desc = 'Add the header of 1st file to 2nd file.')
pAddHeader.input          = "infile1:file, infile2:file"
pAddHeader.output         = "outfile:file:{{in.infile2 | bn}}"
pAddHeader.args.n         = 1
pAddHeader.script         = "file:scripts/common/pAddHeader.bash"


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
pMergeFiles               = Proc(desc = 'Merge files.')
pMergeFiles.input         = "infiles:files"
pMergeFiles.output        = "outfile:file:{{in.infiles[0] | fn}}_etc{{in.infiles[0] | ext}}"
pMergeFiles.args.skip     = 0
pMergeFiles.args.comment  = '#'
pMergeFiles.args.header   = False
pMergeFiles.lang          = params.python.value
pMergeFiles.script        = "file:scripts/common/pMergeFiles.py"

pSplitRows                = Proc(desc = 'Split a file by rows.')
pSplitRows.input          = 'infile:file'
pSplitRows.output         = 'outdir:dir:{{in.infile | bn}}-rows'
pSplitRows.args.skip      = 0
pSplitRows.args.cnames    = True
pSplitRows.args.n         = 8
pSplitRows.lang           = params.python.value
pSplitRows.script         = "file:scripts/common/pSplitRows.py"

