from diot import Diot
from bioprocs.utils import shell2 as shell, logger

infile  = {{i.infile | quote}}
sqlfile = {{i.sqlfile | quote}}
outfile = {{o.outfile | quote}}
argssql = {{args.sql | quote}}
inopts  = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}

if sqlfile:
	with open(sqlfile) as f:
		sql = ' '.join(f.readlines()).strip()
	if argssql:
		logger.warning('`args.sql` is ignored, as `i.sqlfile` is provided.')
else:
	sql = argssql

if not sql:
	raise ValueError('One of `i.sqlfile` and `args.sql` is requied.')

params = {
	'H': inopts.cnames,
	'd': inopts.delimit,
	'e': inopts.encoding,
	'z': (inopts.gz == 'auto' and infile.endswith('.gz')) or inopts.gz is True,
	'D': outopts.delimit if outopts.delimit is not None else inopts.delimit,
	'O': outopts.cnames if outopts.cnames is not None else inopts.cnames,
	'E': outopts.encoding if outopts.encoding is not None else inopts.encoding,
	' ': sql
}
shell.cat(infile).p | shell.q(**params) ^ outfile
