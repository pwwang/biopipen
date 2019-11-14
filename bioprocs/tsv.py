"""TSV file operations"""

from pyppl import Proc, Box
from . import params, rimport
from .utils import fs2name
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pMatrixR():
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
	return Box(
		desc   = 'Operate a matrix and save the new matrix to file',
		lang   = params.Rscript.value,
		input  = "infile:file",
		output = "outfile:file:{{i.infile | bn}}",
		args   = Box(inopts = Box(cnames = True, rnames = True, delimit = "\t", skip = 0),
			params = Box({
				"check.names": "FALSE",
				"quote"      : ""
			}),
			code = []))

@procfactory
def _pTranspose():
	"""
	@name:
		pTranspose
	@description:
		Transpose a matrix
	@input:
		`infile:file`: The input matrix file
	@output:
		`outfile:file`: The Transposed matrix file. Default: `{{i.infile | bn}}`
	@args:
		`inopts`: Input options for input file.
	"""
	pTranspose              = Proc(desc = 'Transpose a matrix')
	pTranspose.input        = 'infile:file'
	pTranspose.output       = 'outfile:file:{{i.infile | bn}}'
	pTranspose.args.inopts  = Box(cnames = True, rnames = True)
	pTranspose.envs.rimport = rimport
	pTranspose.lang         = params.Rscript.value
	pTranspose.script       = "file:scripts/tsv/pTranspose.r"
	return pTranspose

@procfactory
def _pPaired():
	"""
	@name:
		pPaired
	@description:
		Subset each input file and make sure they have paired columns.
	@input:
		`infile1:file`: The first file
		`infile2:file`: The second file
	@outfile:
		`outfile1:file`: The paired file for infile1. Default: `{{i.infile1 | fn}}.paired{{i.infile1 | ext}}`
		`outfile2:file`: The paired file for infile2. Default: `{{i.infile2 | fn}}.paired{{i.infile2 | ext}}`
	@args:
		`inopts1`: reading options for input file1
		`inopts2`: reading options for input file2
	"""
	pPaired        = Proc(desc = "Subset each input file and make sure they have paired columns.")
	pPaired.input  = 'infile1:file, infile2:file'
	pPaired.output = [
		'outfile1:file:{{i.infile1 | fn}}.paired{{i.infile1 | ext}}',
		'outfile2:file:{{i.infile2 | fn}}.paired{{i.infile2 | ext}}'
	]
	pPaired.args.inopts1 = Box(head = True, headCallback = None)
	pPaired.args.inopts2 = Box(head = True, headCallback = None)
	pPaired.lang         = params.python.value
	pPaired.script       = "file:scripts/tsv/pPaired.py"
	return pPaired


@procfactory
def _pCbind():
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
	return pCbind

@procfactory
def _pRbind():
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
	return pRbind

@procfactory
def _pCsplit():
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
	return pCsplit

@procfactory
def _pRsplit():
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
	return pRsplit

@procfactory
def _pTsv():
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
		`outopts`: The output options for outfile
		`row`: A row function to transform/filter the row. Argument is an instance of `TsvRecord`
		`helper`: A helper function for `args.ops`
	"""
	pTsv              = Proc(desc = 'Read, Transform, filter a TSV file.')
	pTsv.input        = "infile:file"
	pTsv.output       = "outfile:file:{{i.infile | fn}}.tsv"
	pTsv.lang         = params.python.value
	pTsv.args.helper  = ''
	pTsv.args.row     = None
	pTsv.args.inopts  = Box(delimit = '\t', comment = '#', skip = 0, cnames = True)
	pTsv.args.outopts = Box(delimit = '\t', cnames = True)
	pTsv.script       = "file:scripts/tsv/pTsv.py"
	return pTsv

@procfactory
def _pTsvColFilter(alias = 'pTsvColSelect'):
	"""
	@input:
		infile : The input file
		colfile: The file with columns, one per line, or a list of columns separated by comma.
			- If this is provided, `args.cols` will be ignored.
	@output:
		outfile: The output file
	@args:
		inopts: The options for reading input file. Default: `Box(cnames = True)`
		keep  : Whether to keep in `args.cols` or to discard
		cols  : The columns used to filter. Could be names or indices(0-based) or a file containing the column names, one per line.
	"""
	return Box(
		desc   = 'Filter a tsv file by columns',
		lang   = params.python.value,
		input  = 'infile:file, colfile:var',
		output = 'outfile:file:{{i.infile | bn}}',
		args   = Box(
			inopts = Box(cnames = True),
			keep = True,
			cols = None,
		)
	)

@procfactory
def _pTsvAggregate():
	"""
	@name:
		pTsvAggregate
	@description:
		Aggregate on columns with a set of records
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.aggr.txt`
			- With columns `args.on` and aggregated results from `args.aggrs`
			- If `args.on` is a function, then the calculated term will be add to the 1st column.
	@args:
		`inopts`: The options to read the input file, Default: `Box(cnames = True)`
		`on`: Aggregate according to which column, Default: `0`
			- It also could column name if `args.inopts = True`
			- The input file has to sorted by this column
			- Or a string of (lambda) function to calculate the term to aggregate on.
		`helper`: Raw codes to give some help for `args.aggrs`
		`aggrs`: The aggregation methods. Required.
			- It's a `Box` with the keys for aggregated results
			- If `args.inopts.cnames = True` then these keys will be output as column names, otherwise ignored
			- You can also combine the aggregation results.
				- For example: `{"sum,mean": lambda rs: [sum(r.value for r in rs), sum(r.value for r in rs)/float(len(rs))]}`
			- We have some built-in aggregation methods:
				- `args.aggrs.Sum = "$sum:3"`          : Get the sum of the 4th column
				- `args.aggrs.Mean = "$mean:Height`    : Get the mean of column "Height"
				- `args.aggrs.Median = "$median:1"`    : Get the median of the 2nd column
				- `args.aggrs.Min = "$min:1"`          : Get the min of the 2nd column
				- `args.aggrs.Max = "$max:1"`          : Get the max of the 2nd column
				- `args.aggrs.Max2 = "$max:2"`         : Get the max of the 3rd column
				- `args.aggrs["Max,Max2"] = "$max:1,2"`: Get the max of the 2nd and 3rd column, respectively
				- `args.aggrs.CombinedP = "$fisher:1"` : Get the combined pvalues for 1st column using fisher'method (`scipy.stats.combine_pvalues`)
	"""
	pTsvAggregate             = Proc(desc = 'Aggregate on columns with a set of records')
	pTsvAggregate.input       = 'infile:file'
	pTsvAggregate.output      = 'outfile:file:{{i.infile | fn2}}.aggr.txt'
	pTsvAggregate.args.inopts = Box(cnames = True)
	pTsvAggregate.args.on     = 0 # which column
	pTsvAggregate.args.aggrs  = Box()
	pTsvAggregate.args.helper = ''
	pTsvAggregate.lang        = params.python.value
	pTsvAggregate.script      = "file:scripts/tsv/pTsvAggregate.py"
	return pTsvAggregate

@procfactory
def _pTsvHeader():
	"""
	@name:
		pTsvHeader
	@description:
		Get the header of a TSV file
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.header.txt`
	@args:
		`inopts`: The options to read input file. Default: `Box(cnames = True)`
		`filter`: The filter for the header. Default: `None`
			- `None`: no filter
			- `lambda cnames: ...` A callback to manipulate colnames.
	"""
	pTsvHeader             = Proc(desc = 'Get the header of a tsv file.')
	pTsvHeader.input       = 'infile:file'
	pTsvHeader.output      = 'outfile:file:{{i.infile | fn2}}.header.txt'
	pTsvHeader.args.inopts = Box(cnames = True)
	pTsvHeader.args.filter = None
	pTsvHeader.lang        = params.python.value
	pTsvHeader.script      = "file:scripts/tsv/pTsvHeader.py"
	return pTsvHeader

@procfactory
def _pTsvReplaceHeader():
	"""
	@name:
		pTsvReplaceHeader
	@description:
		Replace the header of a TSV file
	@input:
		`infile:file`: The input file
		`hfile:file`:  The file containing the headers, one per line.
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | bn}}`
	@args:
		`inopts`: The options to read input file, Default: `Box(cnames = True)`
		`cnames`: The column names or callback, Default: `None`
			- `None`: use the header in `i.hfile`
			- `<list/str/file>`: the header to use if `i.hfile` is not provided
			- `lambda cnames: ...`: The callback to modify header in `i.hfile` if provided, otherwise modify the original header.
	"""
	pTsvReplaceHeader             = Proc(desc = "Replace the header of a tsv file.")
	pTsvReplaceHeader.input       = 'infile:file, hfile:file'
	pTsvReplaceHeader.output      = 'outfile:file:{{i.infile | bn}}'
	pTsvReplaceHeader.args.inopts = Box(cnames = True)
	pTsvReplaceHeader.args.cnames = None
	pTsvReplaceHeader.lang        = params.python.value
	pTsvReplaceHeader.script      = "file:scripts/tsv/pTsvReplaceHeader.py"
	return pTsvReplaceHeader

@procfactory
def _pTsvJoin():
	"""
	@name:
		pTsvJoin
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
		infiles:files: The input files
	@output:
		outfile:file: The output file
	@args:
		inopts: The input options for infile, each can be a list if different for `infiles`:
			- skip   : First N lines to skip. Default: `0`
			- delimit: The delimit. Default: `\t`
			- comment: The comment line mark. Default: `#`
			- cnames : Whether input file has head. Default: `True`
		outopts: The output options:
			- delimit: The delimit. Default: `\t`
			- cnames : Whether to output the head? Default: `False`
			- headCallback: a string callback to organize headers.
				- E.g: `lambda cnames: ','.join(cnames)`
				- The argument `cnames` will be the joint cnames from all readers
				- Could be `True` to use `\t` to assemble the headers
		debug: Save debug information in stderr file. Default: `False`
		match: The match function.
			- Return -1 if matched, otherwise return index the of record that should be read next.
			- E.g: `lambda r1, r2: -1 if r1.id == r2.id else 0 if r1.id < r2.id else 1`
			- By default, will use the first column to compare
		do: The do function. Global vaiable `fout` is available to write results to output file.
			- E.g: `lambda writer, r1, r2: writer.write(r1)`
		helper: Some helper codes.
	"""
	pTsvJoin              = Proc(desc = 'Read files simultaneously.')
	pTsvJoin.input        = 'infiles:files'
	pTsvJoin.output       = 'outfile:file:{{i.infiles[0] | fn}}.etc.joined.txt'
	pTsvJoin.echo         = 0
	pTsvJoin.args.inopts  = Box(delimit = '\t', skip = 0, comment = '#', cnames = True)
	pTsvJoin.args.outopts = Box(delimit = '\t', cnames = False)
	pTsvJoin.args.debug   = False
	pTsvJoin.args.match   = None
	pTsvJoin.args.do      = None
	pTsvJoin.args.helper  = ''
	pTsvJoin.lang         = params.python.value
	pTsvJoin.script       = "file:scripts/tsv/pTsvJoin.py"
	return pTsvJoin

@procfactory
def _pTsvSql():
	"""
	@name:
		pTsvSql
	@description:
		Query tsv file using SQL. (see: http://harelba.github.io/q/examples.html)
	@input:
		`infile:file` : The input tsv file
		`sqlfile:file`: The file containing the SQLs. If provided, `args.sql` will be ignored.
	@output:
		`outfile:file`: The output file
	@args:
		`sql`: If SQL to execute. Use `-` for table name
		`inopts`: Options for input file.
			- `cnames`: Input file has header? Default: `True`
			- `delimit`: The delimit of input file. Default: `\t`
			- `encoding`: Encoding of input file. Default: `UTF-8`
			- `gz`: Whether input file is gzipped. Default: `auto` (detected from file extension)
		`outopts`: Output options.
			- `cnames`: Inherited from `args.inopts`
			- `delimit`: Inherited from `args.inopts`
			- `encoding`: Inherited from `args.inopts`
	@requires:
		[`q`](http://harelba.github.io/q/index.html)
		This process is built on `q 1.7.1`
	"""
	pTsvSql             = Proc(desc = 'Query tsv file using SQL')
	pTsvSql.input       = 'infile:file, sqlfile:file'
	pTsvSql.output      = 'outfile:file:{{i.infile | fn}}.bysql{{i.infile | ext}}'
	pTsvSql.args.sql    = ''
	pTsvSql.args.inopts = Box(
		cnames = True, delimit = "\t", encoding = 'UTF-8', gz = 'auto')
	pTsvSql.args.outopts = Box(
		cnames = None, delimit = None, encoding = None)
	pTsvSql.lang   = params.python.value
	pTsvSql.script = "file:scripts/tsv/pTsvSql.py"
	return pTsvSql

@procfactory
def _pTsvSample():
	"""
	@name:
		pTsvSample
	@description:
		Sample records from a TSV file
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file, Default: `{{i.infile | fn2}}.sampled.txt`
	@args:
		`inopts`   : input options, only skip available, Default: `Box()`
		`n`        : how many records to sample, Default: `10`
		`arsample` : sample program by Alex Reynolds, Default: `<params.arsample>`
		`replace`  : Whether sample with replacement or not, Default: `False`
		`keeporder`: Keep the order of the sampled records as it's in input file, Default: `False`
		`seed`: The seed, Default: `0`
		`params`: Other params for arsample, Default: `Box()`
	"""
	pTsvSample                = Proc(desc = 'Sample records from a TSV file.')
	pTsvSample.input          = 'infile:file'
	pTsvSample.output         = 'outfile:file:{{i.infile | fn2}}.sampled.txt'
	pTsvSample.args.inopts    = Box()
	pTsvSample.args.n         = 10
	pTsvSample.args.arsample  = params.arsample.value
	pTsvSample.args.replace   = False
	pTsvSample.args.keeporder = False
	pTsvSample.args.seed      = 0
	pTsvSample.args.params    = Box()
	pTsvSample.lang           = params.python.value
	pTsvSample.script         = "file:scripts/tsv/pTsvSample.py"
	return pTsvSample

@procfactory
def _pTsvMerge():
	"""
	@name:
		pTsvMerge
	@description:
		Merge files in the input directory
	@input:
		`indir:file`: The input directory
	@output:
		`outfile:file`: The output file
	@args:
		`inopts`: The options for input file. Default: `Box(skip = 0, comment = '#', delimit = '\t')`
		`outopts`: The options for output file. Default: `Box()`
	"""
	pTsvMerge               = Proc(desc = 'Merge files by rows.')
	pTsvMerge.input         = "infiles:files"
	pTsvMerge.output        = "outfile:file:{{i.infiles | fs2name}}"
	pTsvMerge.args.inopts   = Box(skip = 0, comment = '#', delimit = '\t')
	pTsvMerge.args.outopts  = Box()
	# IOError: [Errno 24] Too many open files
	pTsvMerge.args.maxopen  = 100
	pTsvMerge.envs.fs2name  = fs2name
	pTsvMerge.lang          = params.python.value
	pTsvMerge.script        = "file:scripts/tsv/pTsvMerge.py"
	return pTsvMerge

@procfactory
def _pMergeRows():
	"""
	@name:
		pMergeRows
	@description:
		Merge repeated rows
	@input:
		`infile:file`: The input file (has to be sorted by the repeated columns)
	@output:
		`outfile:file`: The output file. Default: `{{i.infile | bn}}`
	@args:
		`inopts`: The options for input file.
			- defaults: skip: 0, comment: #, delimit '\\t'
		`outopts`: The options for output file. Defaults:
			- head: False (not output head line)
			- headPrefix: `#` (The prefix for head line)
			- headDelimit: `\\t` (The delimit for head line)
			- headTransform: `None` (The callback for head line)
			- delimit: `\\t` (The delimit for output line)
		`match`: The function to return a value to decide whether the row is repeated, argument is a `TsvRecord`.
		`do`   : The merge function in python, argument is a list of `TsvRecord`s or `list`s if `args.inopts.ftype` is `nometa`
	"""
	pMergeRows              = Proc(desc = 'Merge repeated rows.')
	pMergeRows.input        = 'infile:file'
	pMergeRows.output       = 'outfile:file:{{i.infile | bn}}'
	pMergeRows.args.inopts  = Box(skip = 0, comment = '#', delimit = '\t')
	pMergeRows.args.outopts = Box(head = False, headPrefix = '', headDelimit = '\t', headTransform = None, delimit = '\t')
	pMergeRows.args.match   = None
	pMergeRows.args.do      = None
	pMergeRows.lang         = params.python.value
	pMergeRows.script       = "file:scripts/tsv/pMergeRows.py"
	return pMergeRows

