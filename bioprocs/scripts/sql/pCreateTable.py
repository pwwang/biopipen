# sqlite:file=/path/to/sqlite.db;...
{{dsnparse}}
{{schemaparse}}

dsn = dsnparse({{i.dsn | quote}})
if not hasattr(dsn, 'scheme'):
	raise Exception('Cannot determine the scheme from DSN string.')

name, fields = schemaparse({{i.schema | quote}}, type = {{args.intype | quote}}, scheme = dsn.scheme, delimit = {{args.delimit | quote}})

if dsn.scheme in ['sqlite', 'sqlite3']:
	assert hasattr(dsn, 'file') and dsn.file != ':memory:'
	import sqlite3
	conn = sqlite3.connect(dsn.file)
	{# Drop the existing table first #}
	{% if args.drop %}
	dropsql = 'DROP TABLE IF EXISTS "%s"' % name
	conn.execute(dropsql)
	{% endif %}

	createsql  = 'CREATE TABLE IF NOT EXISTS "%s" (\n' % name
	for i, field in enumerate(fields):
		fname, ftype, fstat = field
		if i == len(fields) - 1:
			createsql += '  "%s" %s %s\n' % (fname, ftype, fstat)
		else:
			createsql += '  "%s" %s %s,\n' % (fname, ftype, fstat)
	createsql += ') %s' % dsn.after

	stderr.write('Executing SQL:\n')
	stderr.write(createsql + '\n')
	conn.execute(createsql)
	conn.commit()
	conn.close()
else:
	raise NotImplementedError('Only support sqlite3 currently.')


