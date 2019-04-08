import math
from pyppl import Box
from functools import partial
from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
inopts  = {{ args.inopts | repr}}
on      = {{ args.on if isinstance(args.on, int) or args.on.startswith('lambda') else quote(args.on) }}
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
	helper = [helper]

# built-in aggregators
def aggr_sum(rs, idx):
	ret = [0] * len(idx)
	for r in rs:
		for i, ix in enumerate(idx):
			ret[i] += float(r[ix])
	return ret[0] if len(idx) == 1 else ret

def aggr_mean(rs, idx):
	lenrs = float(len(rs))
	sums  = aggr_sum(rs, idx)
	if len(idx) == 1:
		return sums/lenrs
	return [s/lenrs for s in sums]

def get_median(ls):
	ls = list(sorted(float(l) for l in ls))
	lenls = len(ls)
	if lenls % 2 == 1:
		return ls[int((lenls-1)/2.0)]
	return (ls[int((lenls-2)/2.0)] + ls[int(lenls/2.0)])/2.0

def aggr_median(rs, idx):
	if len(idx) == 1:
		return get_median(r[idx[0]] for r in rs)
	lss = [[] for _ in idx]
	for r in rs:
		for i, ix in enumerate(idx):
			lss[i].append(r[ix])
	return [get_median(ls) for ls in lss]

def aggr_min(rs, idx):
	if len(idx) == 1:
		return min(float(r[idx[0]]) for r in rs)
	ret = [None] * len(idx)
	for r in rs:
		for i, ix in enumerate(idx):
			ret[i] = r[ix] if ret[i] is None else min(ret[i], float(r[ix]))
	return ret

def aggr_max(rs, idx):
	if len(idx) == 1:
		return max(float(r[idx[0]]) for r in rs)
	ret = [None] * len(idx)
	for r in rs:
		for i, ix in enumerate(idx):
			ret[i] = r[ix] if ret[i] is None else max(ret[i], float(r[ix]))
	return ret

def aggr_fisher(rs, idx):
	from scipy.stats import combine_pvalues
	if (len(idx) == 1):
		return combine_pvalues([float(r[idx[0]]) for r in rs])[1]
	ret = [None] * len(idx)
	for i, ix in enumerate(idx):
		ret[i] = combine_pvalues([float(r[ix]) for r in rs])[1]
	return ret

builtin = {
	"sum"   : aggr_sum,
	"mean"  : aggr_mean,
	"median": aggr_median,
	"min"   : aggr_min,
	"max"   : aggr_max,
	"fisher": aggr_fisher,
}

helper = [line for line in helper if line]
exec('\n'.join(helper), globals())
aggrs  = Box()
naggrs = {}
{% for col, func in args.aggrs.items() %}
col  = {{col | quote}}
aggrs [col] = {{func if not func.startswith('$') else quote(func)}}
if not callable(aggrs[col]) and not aggrs[col].startswith('$'):
	raise ValueError("I don't know how for {!r}.".format(col))
if not callable(aggrs[col]):
	fn, args = aggrs[col].split(':', 1)
	fn = fn.strip()[1:]
	if fn not in builtin:
		raise ValueError('Unknown builtin aggregation function, expect sum, mean, median, min or max.')
	args = [int(arg) if arg.isdigit() else arg for arg in args.strip().split(',')]
	aggrs[col] = partial(builtin[fn], idx = args)
naggrs[col] = len(alwaysList(col))
{% endfor %}


reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, delimit = inopts.get('delimit', "\t"))
if reader.cnames:
	if isinstance(on, int):
		writer.cnames = [reader.cnames[on]] + alwaysList(','.join(aggrs.keys()))
	elif callable(on):
		writer.cnames = ['AggrGroup'] + alwaysList(','.join(aggrs.keys()))
	elif on not in reader.cnames:
		raise ValueError('{!r} is not a valid column name!'.format(on))
	else:
		writer.cnames = [on] + alwaysList(','.join(aggrs.keys()))
		on = reader.cnames.index(on)
elif callable(on):
	writer.cnames = ['AggrGroup'] + alwaysList(','.join(aggrs.keys()))
elif not isinstance(on, int):
	raise ValueError('The input file does not have column name (args.inopts.cnames = False?).')
else:
	writer.cnames = ['ROWNAME'] + alwaysList(','.join(aggrs.keys()))

writer.writeHead()

refcol = None
rows   = []
for row in reader:
	if callable(on):
		col = on(row)
	else:
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



