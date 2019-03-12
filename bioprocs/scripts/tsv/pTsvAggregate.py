from pyppl import Box
from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
inopts  = {{ args.inopts | repr}}
on      = {{ args.on | repr}}
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
	helper = [helper]

helper = [line for line in helper if line]
exec('\n'.join(helper), globals())
aggrs  = Box()
naggrs = {}
{% for col, func in args.aggrs.items() %}
col  = {{col | quote}}
aggrs [col] = {{func}}
if not callable(aggrs[col]):
	raise ValueError("I don't know how for {!r}.".format(col))
naggrs[col] = len(alwaysList(col))
{% endfor %}


reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, delimit = inopts.get('delimit', "\t"))
if reader.cnames:
	if isinstance(on, int):
		writer.cnames = [reader.cnames[on]] + alwaysList(','.join(aggrs.keys()))
	elif on not in reader.cnames:
		raise ValueError('{!r} is not a valid column name!'.format(on))
	else:
		writer.cnames = [on] + alwaysList(','.join(aggrs.keys()))
		on = reader.cnames.index(on)
elif not isinstance(on, int):
	raise ValueError('The input file does not have column name (args.inopts.cnames = False?).')
else:
	writer.cnames = ['ROWNAME'] + alwaysList(','.join(aggrs.keys()))

writer.writeHead()

refcol = None
rows   = []
for row in reader:
	col = row[on]
	if refcol is None:
		refcol = col
		rows.append(row)
	elif col == refcol:
		rows.append(row)
	else:
		out = [refcol]
		for c, a in aggrs.items():
			if naggrs[c] > 1:
				out.extend(a(rows))
			else:
				out.append(a(rows))
		writer.write(out)
		refcol = col
		rows   = [row]
out = [refcol]
for c, a in aggrs.items():
	if naggrs[c] > 1:
		out.extend(a(rows))
	else:
		out.append(a(rows))
writer.write(out)
writer.close()



