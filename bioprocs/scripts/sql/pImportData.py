{{dsnparse}}
{{schemaparse}}

from sys import stderr
dsn = dsnparse({{i.dsn | quote}})
if not hasattr(dsn, 'scheme'):
	raise Exception('Cannot determine the scheme from DSN string.')

# maybe needed for other scheme
name, fields = schemaparse({{i.datafile | quote}}, type = 'data', scheme = dsn.scheme, delimit = {{args.delimit | quote}})
name = {{i.datafile | tablename | quote}}

if dsn.scheme in ['sqlite', 'sqlite3']:
	assert hasattr(dsn, 'file') and dsn.file != ':memory:'
	import sqlite3
	conn = sqlite3.connect(dsn.file)
	conn.text_factory = str

	# drop table
	{% if args.drop %}
	dropsql = 'DROP TABLE IF EXISTS "%s"' % name
	conn.execute(dropsql)
	conn.commit()
	{% endif %}

	# create table
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

	insertsql = 'INSERT INTO "%s" ("%s") VALUES (%s);' % (name, '", "'.join([
		field[0] for field in fields
	]), ', '.join(['?' for _ in fields]))
	data = []
	with open({{i.datafile | quote}}, 'r') as f:
		f.readline() # skip header
		for line in f:
			line = line.strip('\n')
			if not line: continue
			line = line.split({{args.delimit | quote}})
			data.append(tuple(line))
	if data:
		stderr.write('Executing SQL:\n')
		stderr.write(insertsql + '\n')
		conn.executemany(insertsql, data)
	conn.commit()
	conn.close()
else:
	raise NotImplementedError('Only support sqlite3 currently.')