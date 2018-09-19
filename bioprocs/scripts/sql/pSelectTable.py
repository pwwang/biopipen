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
	cur  = conn.cursor()
	cur.execute({{args.sql | quote}})
	with open({{out.outfile | quote}}, 'w') as fout:
		keys = [d[0] for d in cur.description]
		fout.write('\t'.join(keys) + '\n')
		for row in cur:
			outs = [str(r) for r in row]
			fout.write('\t'.join(outs) + '\n')
	cur.close()
	conn.close()
else:
	raise NotImplementedError('Only support sqlite3 currently.')