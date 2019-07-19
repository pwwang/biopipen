"""
Some sql database utilities
"""

from pyppl import Proc

from . import params
from .utils import sql
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pCreateTable():
	"""
	@name:
		pCreateTable
	@description:
		Create tables in the database
	@input:
		`dsn`: The dsn to connect to the database
			- currently support `sqlite:file=...`
		`schema:file`: The schema file
			- could be a pure schema file:
			```
			Field	Type	Statement
			ID	INT	PRIMARY KEY
			...
			```
			- or a data file with header
	@output:
		`dsn`: The dsn
	@args:
		`intype`: The input file schema file or a data file. Default: `schema`
		`drop`:  Force creating the table (drop the pre-existing table)
		`delimit`:The delimit of input file. Default: `\\t`
	"""
	pCreateTable                  = Proc(desc = 'Create tables in the database.')
	pCreateTable.input            = "dsn, schema:file"
	pCreateTable.output           = "dsn:var:{{i.dsn}}"
	pCreateTable.args.intype      = 'schema'
	pCreateTable.args.drop        = False
	pCreateTable.args.delimit     = "\t"
	pCreateTable.envs.dsnparse    = sql.dsnparse
	pCreateTable.envs.schemaparse = sql.schemaparse
	pCreateTable.lang             = params.python.value
	pCreateTable.script           = 'file:scripts/sql/pCreateTable.py'
	return pCreateTable

@procfactory
def _pImportData():
	"""
	@name:
		pImportData
	@description:
		Create tables and import the data
	@input:
		`dsn`: The dsn to connect to the database
			- currently support `sqlite:file=...`
		`datafile:file`: The schema file
			- must have header
	@output:
		`dsn`: The dsn
	@args:
		`delimit`:The delimit of input file. Default: `\\t`
	"""
	pImportData                  = Proc(desc = 'Create tables and import the data')
	pImportData.input            = "dsn, datafile:file"
	pImportData.output           = "dsn:var:{{i.dsn}}"
	pImportData.args.delimit     = "\t"
	pImportData.args.drop        = False
	pImportData.envs.dsnparse    = sql.dsnparse
	pImportData.envs.schemaparse = sql.schemaparse
	pImportData.envs.tablename   = lambda fn: __import__('os').path.basename(fn).split('.')[0]
	pImportData.lang             = params.python.value
	pImportData.script           = 'file:scripts/sql/pImportData.py'
	return pImportData

@procfactory
def _pUpdateTable():
	"""
	@name:
		pUpdateTable
	@description:
		Update table using sql.
	@input:
		`dsn`: The dsn to connect to the database
			- currently support `sqlite:file=...`
	@output:
		`dsn`: The dsn
	@args:
		`sql`: The sql to update the table (list)
	"""
	pUpdateTable               = Proc(desc = 'Update table using sql.')
	pUpdateTable.input         = 'dsn'
	pUpdateTable.output        = 'dsn:var:{{i.dsn}}'
	pUpdateTable.args.sql      = []
	pUpdateTable.envs.dsnparse = sql.dsnparse
	pUpdateTable.lang          = params.python.value
	pUpdateTable.script        = 'file:scripts/sql/pUpdateTable.py'
	return pUpdateTable

@procfactory
def _pSelectTable():
	"""
	@name:
		pSelectTable
	@description:
		Select data from table and dump it.
	@input:
		`dsn`: The dsn to connect to the database
			- currently support `sqlite:file=...`
	@output:
		`outfile:file`: The dumped file
	@args:
		`sql`: The sql to select data from the table (list)
	"""
	pSelectTable               = Proc(desc = 'Select data from table and dump it.')
	pSelectTable.input         = 'dsn'
	pSelectTable.output        = 'outfile:file:sqlSelectTable-{{i.dsn, args.sql | lambda x,y: __import__("hashlib").md5(str(x) + str(y)).hexdigest()[:16]}}.dumped.txt'
	pSelectTable.args.sql      = []
	pSelectTable.envs.dsnparse = sql.dsnparse
	pSelectTable.lang          = params.python.value
	pSelectTable.script        = 'file:scripts/sql/pSelectTable.py'
	return pSelectTable

