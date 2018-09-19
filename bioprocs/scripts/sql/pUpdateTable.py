{{dsnparse}}

from sys import stderr
dsn = dsnparse({{i.dsn | quote}})
if not hasattr(dsn, 'scheme'):
	raise Exception('Cannot determine the scheme from DSN string.')

if dsn.scheme in ['sqlite', 'sqlite3']:
	assert hasattr(dsn, 'file') and dsn.file != ':memory:'
	import sqlite3
	conn = sqlite3.connect(dsn.file)
	conn.text_factory = str
	for sql in {{args.sql}}:
		conn.execute(sql)
		stderr.write('Executing SQL:\n')
		stderr.write(createsql + '\n')
	conn.commit()
	conn.close()
else:
	raise NotImplementedError('Only support sqlite3 currently.')