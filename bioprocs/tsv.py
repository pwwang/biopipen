from pyppl import Proc, Box
from . import params

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
pRbind.args.cnames    = True # or [True, True, False] corresponding to the file order
pRbind.args.rnames    = True
pRbind.args.na        = 'NA'
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
	`inmeta`: The meta data for input file. Could be:
		- A list of columns names,
		- A list of OrderedDict items with column name and its description
		- An OrderedDict (no dict (won't be checked!!), because key order won't keep)
		- A string for predefined reader classes. For example, bed, bed12, bedpe, bedx, etc.
	`outmeta`: The meta data for output file. Same as inmeta, but you can specify part of them to just keep some columns. You may also use a different column name, but you have to specify values for each row in `args.ops` function
	`ops`: A ops function to transform the row. Argument is an instance of `readRecord`
	`opshelper`: A helper function for `args.ops`
	`inopts`: The options for reader. Default: `Box(delimit = '\\t', comment = '#', skip = 0)`
	`outdem`: The output delimiter. It won't work if it is a predefined write class
	`incom`: The prefix of comments. It won't work if it is a predefined read class
	`outdem`: The prefix of comments. It won't work if it is a predefined write class
	`omprefix`: The prefix for output metadata. Defualt: '##META/'
	`hdprefix`: The prefix for output header. Default: '#'
	`writemeta`: Whether report meta to the output file. Default: True
	`writehead`: Whether report header to the output file. Default: True
"""
pTsv                = Proc(desc = 'Read, Transform, filter a TSV file.')
pTsv.input          = "infile:file"
pTsv.output         = "outfile:file:{{in.infile | fn}}.tsv"
pTsv.lang           = params.python.value
pTsv.args.inmeta    = None
pTsv.args.outmeta   = None
pTsv.args.opshelper = ''
pTsv.args.ops       = None
pTsv.args.inopts    = Box(delimit = '\t', comment = '#', skip = 0)
pTsv.args.outopts   = Box(delimit = '\t', metaprefix = '##META/', headprefix = '#', meta = True, head = True)
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
pSimRead.args.skip    = 0
pSimRead.args.usehead = None
pSimRead.args.delimit = '\t'
pSimRead.args.gzip    = 'auto'
pSimRead.args.match   = None
pSimRead.args.do      = None
pSimRead.args.data    = {}
pSimRead.lang         = params.python.value
pSimRead.script       = "file:scripts/tsv/pSimRead.py"