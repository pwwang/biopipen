from pyppl import Proc, Box
from . import params
from .utils import helpers, txt

"""
@name:
	pMatrix
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
pMatrix             = Proc(desc = 'Operate a matrix and save the new matrix to file.')
pMatrix.input       = "infile:file"
pMatrix.output      = "outfile:file:{{in.infile | bn}}"
pMatrix.args.cnames = True
pMatrix.args.rnames = True
pMatrix.args.code   = []
pMatrix.lang        = params.Rscript.value
pMatrix.script      = "file:scripts/matrix/pMatrix.r"

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
pCbind.output         = 'outfile:file:{{in.infiles[0] | bn}}'
pCbind.args.cnames    = True # or [True, True, False] corresponding to the file order
pCbind.args.rnames    = True
pCbind.args.na        = 'NA'
pCbind.envs.cbindfill = helpers.cbindfill.r
pCbind.lang           = params.Rscript.value
pCbind.script         = "file:scripts/matrix/pCbind.r"

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
pRbind.args.cnames    = True # or [True, True, False] corresponding to the file order
pRbind.args.rnames    = True
pRbind.args.na        = 'NA'
pRbind.envs.rbindfill = helpers.rbindfill.r
pRbind.lang           = params.Rscript.value
pRbind.script         = "file:scripts/matrix/pRbind.r"

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
pCsplit.output      = 'outdir:dir:{{in.infile | fn}}.splits'
pCsplit.args.cnames = True
pCsplit.args.rnames = True
pCsplit.lang        = params.Rscript.value
pCsplit.script      = "file:scripts/matrix/pCsplit.r"

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
pRsplit.output        = 'outdir:dir:{{in.infile | fn}}.splits'
pRsplit.args.cnames   = True
pRsplit.args.rnames   = True
pRsplit.lang          = params.Rscript.value
pRsplit.script        = "file:scripts/matrix/pRsplit.r"

"""
@name:
	pTxtFilter
@description:
	Filter a tab-delimit file (txt/tsv file)
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`cols`:    The cols to remain. Could be either index or column name. Default: `[] (all columns)`
	`rfilter`: The row filter. Default: `None` (don`t filter)
		- Should be a string of a lambda function with parameter of a list of fields.
		- Note that the fields should be the ones remained after column filtering.
		- Rows remained when this function returns `True`.
	`header`:  Whether input file has a header. Default: `True`
	`skip`:    Skip first serveral lines. Default: 0
	`delimit`: The delimit. Default: `\t`
"""
pTxtFilter                 = Proc(desc = 'Filter a txt(tsv) file by columns and rows.')
pTxtFilter.input           = "infile:file"
pTxtFilter.output          = "outfile:file:{{in.infile | bn}}"
pTxtFilter.lang            = params.python.value
pTxtFilter.args.cols       = []
pTxtFilter.args.rfilter    = None
pTxtFilter.args.header     = True
pTxtFilter.args.skip       = 0
pTxtFilter.args.data       = {}
pTxtFilter.args.delimit    = "\t"
pTxtFilter.args.outdelimit = "\t"
pTxtFilter.envs.txtFilter  = txt.filter.py
pTxtFilter.script          = "file:scripts/matrix/pTxtFilter.py"

"""
@name:
	pTxtTransform
@description:
	Transform a tab-delimit file (txt/tsv file)
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`cols`:    The cols to remain. Could be either index or column name. Default: `[] (all columns)`
	`transform`: The row transformer. Default: `None` (don`t transform)
		- Should be a string of a lambda function with parameter of a list of fields.
		- Note that the fields should be the ones remained after column filtering.
		- Must return a list of strings(fields).
	`header`:  Whether input file has a header. Default: `True`
	`skip`:    Skip first serveral lines. Default: 0
	`delimit`: The delimit. Default: `\t`
"""
pTxtTransform                   = Proc(desc = 'Transform a txt(tsv) file.')
pTxtTransform.input             = "infile:file"
pTxtTransform.output            = "outfile:file:{{in.infile | bn}}"
pTxtTransform.lang              = params.python.value
pTxtTransform.args.cols         = []
pTxtTransform.args.data         = {}
pTxtTransform.args.transform    = None
pTxtTransform.args.header       = True
pTxtTransform.args.skip         = 0
pTxtTransform.args.delimit      = "\t"
pTxtTransform.args.outdelimit   = "\t"
pTxtTransform.envs.txtTransform = txt.transform.py
pTxtTransform.script            = "file:scripts/matrix/pTxtTransform.py"

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
pSimRead.args.skip    = 0
pSimRead.args.usehead = None
pSimRead.args.delimit = '\t'
pSimRead.args.gzip    = 'auto'
pSimRead.args.match   = None
pSimRead.args.do      = None
pSimRead.args.data    = {}
pSimRead.lang         = params.python.value
pSimRead.script       = "file:scripts/matrix/pSimRead.py"