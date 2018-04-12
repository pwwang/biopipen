from pyppl import Proc, Box
from . import params, rimport

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
	`cnames`: Whether the input file has cnames. Default: True
	`rnames  `: Whether the input file has rnames  . Default: 1
	`code`: The R code to operating the matrix. (the matrix is read in variable `mat`)
"""
pMatrixR             = Proc(desc = 'Operate a matrix and save the new matrix to file.')
pMatrixR.input       = "infile:file"
pMatrixR.output      = "outfile:file:{{in.infile | bn}}"
pMatrixR.args.cnames = True
pMatrixR.args.rnames = True
pMatrixR.args.params = Box({
	"sep"        : "\t",
	"check.names": "FALSE",
	"quote"      : ""
})
pMatrixR.args.code   = []
pMatrixR.lang        = params.Rscript.value
pMatrixR.script      = "file:scripts/tsv/pMatrixR.r"

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
	`cnames`: Whether the input file has cnames. Default: True
		- or [True, True, False] corresponding to the file order
	`rnames  `: Whether the input file has rnames  . Default: 1
	`miss`: Replacement for missing values. Default: `NA`
"""
pCbind                = Proc(desc = 'Cbind the rest of files to the first file.')
pCbind.input          = 'infiles:files'
pCbind.output         = 'outfile:file:{{in.infiles[0] | fn2}}.cbound.txt'
pCbind.args.cnames    = True # or [True, True, False] corresponding to the file order
pCbind.args.rnames    = True
pCbind.args.na        = 'NA'
pCbind.lang           = params.Rscript.value
pCbind.script         = "file:scripts/tsv/pCbind.r"

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
	`cnames`: Whether the input file has cnames. Default: True
		- or [True, True, False] corresponding to the file order
	`rnames  `: Whether the input file has rnames  . Default: 1
	`miss`: Replacement for missing values. Default: `NA`
"""
pRbind                = Proc(desc = 'Rbind the rest of files to the first file.')
pRbind.input          = 'infiles:files'
pRbind.output         = 'outfile:file:{{in.infiles[0] | bn}}'
pRbind.args.inopts    = Box(
	cnames = True, # or [True, True, False] corresponding to the file order
	rnames = True
)
pRbind.args.na        = 'NA'
pRbind.envs.rimport   = rimport
pRbind.lang           = params.Rscript.value
pRbind.script         = "file:scripts/tsv/pRbind.r"

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
	`cnames`: Whether the input file has cnames. Default: True
		- or [True, True, False] corresponding to the file order
	`rnames  `: Whether the input file has rnames  . Default: 1
"""
pCsplit             = Proc(desc = 'Split the columns of input file into different files.')
pCsplit.input       = 'infile:file'
pCsplit.output      = 'outdir:dir:{{in.infile | fn}}.csplits'
pCsplit.args.cnames = True
pCsplit.args.rnames = True
pCsplit.args.n      = 1
pCsplit.lang        = params.Rscript.value
pCsplit.script      = "file:scripts/tsv/pCsplit.r"

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
	`cnames`: Whether the input file has cnames. Default: True
		- or [True, True, False] corresponding to the file order
	`rnames  `: Whether the input file has rnames  . Default: 1
"""
pRsplit               = Proc(desc = 'Rbind the rest of files to the first file.')
pRsplit.input         = 'infile:file'
pRsplit.output        = 'outdir:dir:{{in.infile | fn}}.rsplits'
pRsplit.args.cnames   = True
pRsplit.args.rnames   = True
pRsplit.args.n        = 1
pRsplit.lang          = params.Rscript.value
pRsplit.script        = "file:scripts/tsv/pRsplit.r"

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
		- `delimit`: The delimit. Default: `\\t`
		- `comment`: The comment sign. Default: `#`
		- `skip`: First N lines to skip. Default: `0`
		- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict). If not specified, metadata will be generated automatically.
	`outopts`: The output options for outfile:
		- `delimit`: The delimit for records. Default: `\\t`
		- `head`: Output header or not. Default: `False`
		- `headDelimit`: The delimit for header. Default: `\\t`
		- `headPrefix`: The prefix for header. Default: ``
		- `headTransform`: The transformer for header. Default: `None`
		- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict, '+' as an element or key is allowed to indicate extra meta from the reader). If not specified, metadata will be borrowed from the reader. 
	`ops`: A ops function to transform the row. Argument is an instance of `readRecord`
	`opshelper`: A helper function for `args.ops`
"""
pTsv                = Proc(desc = 'Read, Transform, filter a TSV file.')
pTsv.input          = "infile:file"
pTsv.output         = "outfile:file:{{in.infile | fn}}.tsv"
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
	`skip`: argument skip for each file
	`delimit`: argument delimit for each file
	`usehead`: The header from which input file will be used for output file.
		- Default: None (Don't write header)
	`gzip`: argument gzip for each file
	`match`: The match function. 
	`do`: The do function. Global vaiable `fout` is available to write results to output file.
@requires:
	[`python-simread`](https://github.com/pwwang/simread)
"""
pSimRead              = Proc(desc = 'Read files simultaneously.')
pSimRead.input        = 'infiles:files'
pSimRead.output       = 'outfile:file:{{in.infiles[0] | fn}}.etc.simread.txt'
pSimRead.args.inopts  = Box(delimit = '\t', skip = 0, comment = '#', ftype = 'nometa')
pSimRead.args.outopts = Box(delimit = '\t', headPrefix = '', headDelimit = '\t', headTransform = None, head = False, ftype = '', cnames = [])
pSimRead.args.usemeta = None
pSimRead.args.match   = None
pSimRead.args.do      = None
pSimRead.args.helper  = ''
pSimRead.lang         = params.python.value
pSimRead.script       = "file:scripts/tsv/pSimRead.py"