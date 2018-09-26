from pyppl import Proc, Box
from . import params, rimport
from .utils import fs2name

"""
@name:
	pMatrixR
@description:
	Operate a matrix and save the new matrix to file.
@input:
	`infile:file`: The input file containing the matrix
@output:
	`outfile:file`: The output matrix
@args:
	`inopts`: The input options for infile:
		- `cnames`: Whether the input file has cnames. Default: True
		- `rnames  `: Whether the input file has rnames. Default: True
		- `delimit`: The delimit. Default: `\t`
		- `skip`: First N lines to skip. Default: `0`
	`params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
	`code`: The R code to operating the matrix. (the matrix is read in variable `mat`)
"""
pMatrixR             = Proc(desc = 'Operate a matrix and save the new matrix to file.')
pMatrixR.input       = "infile:file"
pMatrixR.output      = "outfile:file:{{i.infile | bn}}"
pMatrixR.args.inopts = Box(cnames = True, rnames = True, delimit = "\t", skip = 0)
pMatrixR.args.params = Box({
	"check.names": "FALSE",
	"quote"      : ""
})
pMatrixR.args.code    = []
pMatrixR.envs.rimport = rimport
pMatrixR.lang         = params.Rscript.value
pMatrixR.script       = "file:scripts/tsv/pMatrixR.r"

"""
@name:
	pCbind
@description:
	Cbind the rest of files to the first file.
@input:
	`infiles:files`: The input files
@output:
	`outfile:file`: The output matrix
@args:
	`inopts`: The input options for infile:
		- `cnames`: Whether the input file has cnames. Default: True
		   - or [True, True, False] corresponding to the file order
		- `rnames  `: Whether the input file has rnames. Default: True
		- `delimit`: The delimit. Default: `\t`
		- `skip`: First N lines to skip. Default: `0`
	`params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
	`fn2cname`: The function (r) used to convert file name to column name.
	`fill`: Do `cbind.fill` instead of `cbind`. Default: `True`
		- Set it to `False` if the row names are in the same order
	`na`: Replacement for missing values. Default: `NA`
"""
pCbind             = Proc(desc = 'Cbind the rest of files to the first file.')
pCbind.input       = 'infiles:files'
pCbind.output      = 'outfile:file:{{i.infiles | fs2name}}.cbound.txt'
pCbind.args.inopts = Box(
	cnames  = True, # or [True, True, False] corresponding to the file order
	rnames  = True,
	delimit = "\t",
	skip    = 0
)
pCbind.args.params = Box({
	"check.names": "FALSE",
	"quote"      : ""
})
pCbind.args.na       = 'NA'
pCbind.args.fn2cname = 'function(fn) fn' # if 
pCbind.args.fill     = True
pCbind.envs.fs2name  = fs2name
pCbind.envs.rimport  = rimport
pCbind.lang          = params.Rscript.value
pCbind.script        = "file:scripts/tsv/pCbind.r"

"""
@name:
	pRbind
@description:
	Rbind the rest of files to the first file.
@input:
	`infiles:files`: The input files
@output:
	`outfile:file`: The output matrix
@args:
	`inopts`: The input options for infile:
		- `cnames`: Whether the input file has cnames. Default: True
		   - or [True, True, False] corresponding to the file order
		- `rnames  `: Whether the input file has rnames. Default: True
		- `delimit`: The delimit. Default: `\t`
		- `skip`: First N lines to skip. Default: `0`
	`params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
	`na`: Replacement for missing values. Default: `NA`
	`fn2rname`: The function (r) used to convert file name to row name.
	`fill`: Do `rbind.fill` instead of `rbind`. Default: `True`
		- Set it to `False` if the row names are in the same order
	`na`: Replacement for missing values. Default: `NA`
"""
pRbind                = Proc(desc = 'Rbind the rest of files to the first file.')
pRbind.input          = 'infiles:files'
pRbind.output         = 'outfile:file:{{i.infiles[0] | bn}}'
pRbind.args.inopts    = Box(
	cnames  = True, # or [True, True, False] corresponding to the file order
	rnames  = True,
	delimit = "\t",
	skip    = 0
)
pRbind.args.params = Box({
	"check.names": "FALSE",
	"quote"      : ""
})
pRbind.args.na       = 'NA'
pRbind.args.fn2rname = 'function(fn) fn'
pRbind.args.fill     = True
pRbind.envs.fs2name  = fs2name
pRbind.envs.rimport  = rimport
pRbind.lang          = params.Rscript.value
pRbind.script        = "file:scripts/tsv/pRbind.r"

"""
@name:
	pCsplit
@description:
	Split a matrix by columns and save them into files.
@input:
	`infile:file`: The input file
@output:
	`outdir:dir`: The directory containing the output column files
@args:
	`inopts`: The input options for infile:
		- `cnames`: Whether the input file has cnames. Default: True
		- `rnames  `: Whether the input file has rnames. Default: True
		- `delimit`: The delimit. Default: `\t`
		- `skip`: First N lines to skip. Default: `0`
	`params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
	`size`: The chunk size (how many columns to split into one file). Default: `1`
"""
pCsplit             = Proc(desc = 'Split the columns of input file into different files.')
pCsplit.input       = 'infile:file'
pCsplit.output      = 'outdir:dir:{{i.infile | fn}}.csplits'
pCsplit.args.inopts = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t",
	skip    = 0
)
pCsplit.args.params = Box({
	"check.names": "FALSE",
	"quote"      : ""
})
pCsplit.args.size    = 1
pCsplit.envs.rimport = rimport
pCsplit.lang         = params.Rscript.value
pCsplit.script       = "file:scripts/tsv/pCsplit.r"

"""
@name:
	pRsplit
@description:
	Split a matrix by rows and save them into files.
@input:
	`infile:file`: The input file
@output:
	`outdir:dir`: The directory containing the output row files
@args:
	`inopts`: The input options for infile:
		- `cnames`: Whether the input file has cnames. Default: True
		- `rnames  `: Whether the input file has rnames. Default: True
		- `delimit`: The delimit. Default: `\t`
		- `skip`: First N lines to skip. Default: `0`
	`params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`
	`size`: The chunk size (how many rows to split into one file). Default: `1`
"""
pRsplit             = Proc(desc = 'Rbind the rest of files to the first file.')
pRsplit.input       = 'infile:file'
pRsplit.output      = 'outdir:dir:{{i.infile | fn}}.rsplits'
pRsplit.args.inopts = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t",
	skip    = 0
)
pRsplit.args.params = Box({
	"check.names": "FALSE",
	"quote"      : ""
})
pRsplit.args.size    = 1
pRsplit.envs.rimport = rimport
pRsplit.lang         = params.Rscript.value
pRsplit.script       = "file:scripts/tsv/pRsplit.r"

"""
@name:
	pTsv
@description:
	Read, Transform, filter a TSV file.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`inopts`: The input options for infile:
		- `delimit`: The delimit. Default: `\t`
		- `comment`: The comment sign. Default: `#`
		- `skip`: First N lines to skip. Default: `0`
		- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict). If not specified, metadata will be generated automatically.
	`outopts`: The output options for outfile:
		- `delimit`: The delimit for records. Default: `\t`
		- `head`: Output header or not. Default: `False`
		- `headDelimit`: The delimit for header. Default: `\t`
		- `headPrefix`: The prefix for header. Default: ``
		- `headTransform`: The transformer for header. Default: `None`
		- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict, '+' as an element or key is allowed to indicate extra meta from the reader). If not specified, metadata will be borrowed from the reader. 
	`ops`: A ops function to transform the row. Argument is an instance of `readRecord`
	`opshelper`: A helper function for `args.ops`
"""
pTsv                = Proc(desc = 'Read, Transform, filter a TSV file.')
pTsv.input          = "infile:file"
pTsv.output         = "outfile:file:{{i.infile | fn}}.tsv"
pTsv.lang           = params.python.value
pTsv.args.opshelper = ''
pTsv.args.ops       = None
pTsv.args.inopts    = Box(delimit = '\t', comment = '#', skip = 0, ftype = '', cnames = [])
pTsv.args.outopts   = Box(delimit = '\t', headPrefix = '', headDelimit = '\t', headTransform = None, head = False, ftype = '', cnames = [])
pTsv.script         = "file:scripts/tsv/pTsv.py"

"""
@name:
	pSimRead
@description:
	Read files simultaneously.
	NOTE: only one file allows multiple lines with same value to compare, and that file should be the first one. For example: 
	```
	File1:
	1	1
	1	2
	1	3
	File2:
	1	1
	2	2
	3	3
	```
	If you compare the first column, File1 has to put at the begining for input.
@input:
	`infiles:files`: The input files
@output:
	`outfile:file`: The output file
@args:
	`inopts`: The input options for infile:
		- `skip`   : First N lines to skip. Default: `0`
		- `delimit`: The delimit. Default          : `\t`
		- `comment`: The comment line mark. Default: `#`
		- `ftype  `: The type of the file. Default : `nometa`
		- `head`: Whether input file has head, only for `nometa`. Default: `True`
	`outputs`:
		- `delimit`      : The delimit. Default                    : `\t`
		- `headPrefix`   : The prefix for the head line. Default   : ``
		- `headDelimit`  : The delimiter for the head line. Default: `\t`
		- `headTransform`: The transformer for the head line.
		- `head`         : Whether to output the head? Default     : `False`
		- `ftype`        : The type of the output. Default         : `nometa`
		- `cnames`       : The extra column names. Default         : `[]`
	`usemeta`: The header from which input file will be used for output file.
		- Default: None (Don't write header)
		- 1-based.
	`match`: The match function. 
	`do`: The do function. Global vaiable `fout` is available to write results to output file.
	`helper`: Some helper codes.
@requires:
	[`python-simread`](https://github.com/pwwang/simread)
"""
pSimRead              = Proc(desc = 'Read files simultaneously.')
pSimRead.input        = 'infiles:files'
pSimRead.output       = 'outfile:file:{{i.infiles[0] | fn}}.etc.simread.txt'
pSimRead.args.inopts  = Box(delimit = '\t', skip = 0, comment = '#', ftype = 'nometa', head = True)
pSimRead.args.outopts = Box(delimit = '\t', headPrefix = '', headDelimit = '\t', headTransform = None, head = True, ftype = 'nometa', cnames = [])
pSimRead.args.usemeta = None
pSimRead.args.match   = None
pSimRead.args.do      = None
pSimRead.args.helper  = ''
pSimRead.lang         = params.python.value
pSimRead.script       = "file:scripts/tsv/pSimRead.py"

"""
@name:
	pTsvJoin
@description:
	Alias of pSimRead
"""
pTsvJoin = pSimRead.copy()

"""
@name:
	pMergeFiles
@description:
	Merge files in the input directory
@input:
	`indir:file`: The input directory
@output:
	`outfile:file`: The output file
@args:
	`inopts`: The options for input file.
		- defaults: skip: 0, comment: #, delimit '\\t'
	`outopts`: The options for output file. Defaults:
		- head: False (not output head line)
		- headPrefix: `#` (The prefix for head line)
		- headDelimit: `\\t` (The delimit for head line)
		- headTransform: `None` (The callback for head line)
		- delimit: `\\t` (The delimit for output line)
"""
pMergeFiles               = Proc(desc = 'Merge files.')
pMergeFiles.input         = "infiles:files"
pMergeFiles.output        = "outfile:file:{{i.infiles | lambda x: x[0] if x else 'nothing' |  fn}}.etc{{i.infiles | lambda x: x[0] if x else '' | ext}}"
pMergeFiles.args.inopts   = Box(skip = 0, comment = '#', delimit = '\t')
pMergeFiles.args.outopts  = Box(head = False, headPrefix = '', headDelimit = '\t', headTransform = None, delimit = '\t')
# IOError: [Errno 24] Too many open files
pMergeFiles.args.maxopen  = 100
pMergeFiles.lang          = params.python.value
pMergeFiles.script        = "file:scripts/tsv/pMergeFiles.py"