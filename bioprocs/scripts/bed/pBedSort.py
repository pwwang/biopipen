from os import path
from pyppl import Box
from bioprocs.utils import shell

infile    = {{ i.infile | quote}}
outfile   = {{ o.outfile | quote }}
joboutdir = {{job.outdir | quote}}
tool      = {{ args.tool | quote}}
bedtools  = {{ args.bedtools | quote}}
bedops    = {{ args.bedops | quote}}
argsmem   = {{ args.mem | quote}}
sortby    = {{ args.by | quote}}
unique    = {{ args.unique | bool}}
params    = {{ args.params | repr}}
tmpdir    = {{ args.tmpdir | quote}}
chrorder  = {{ args.chrorder if not args.chrorder or isinstance(args.chrorder, (tuple, list)) else args.chrorder.split(',') | repr }}

shell.TOOLS.bedtools = bedtools
shell.TOOLS.bedops   = bedtools

def run_sort():
	params.T = tmpdir
	params.S = argsmem
	params.u = unique
	if sortby == 'coord':
		params.k = ['1,1', '2,2n']
	else:
		params.k = '4'
	shell.grep('^#', infile, _stdout = outfile)
	shell.Shell().grep(v = '^#', _ = infile).pipe(dash = '-', equal = ' ', duplistkey = True).sort(__stdout = outfile, **params).run()

def run_bedops():
	params['max-mem'] = args.mem
	params.tmpdir = tmpdir
	params._ = infile
	shell.grep('^#', infile, _stdout = outfile)
	c = shell.Shell(equal = ' ')
	if unique:
		c.bedops(**params).pipe().uniq(__stdout = outfile).run()
	else:
		params.__stdout = outfile
		c.bedops(**params).run()
	
def run_bedtools():
	params.i = infile
	shell.grep('^#', infile, _stdout = outfile)
	c = shell.Shell(dash = '-', equal = ' ', subcmd = True).bedtools
	if unique:
		c.sort(**params).pipe().uniq(__stdout = outfile).run()
	else:
		params.__stdout = outfile
		c.sort(**params).run()

def attach_order():
	if tool != 'sort':
		return
	mappings = dict(zip(chrorder, sorted(chrorder)))
	global infile
	tmpfile = path.join(joboutdir, path.basename(infile) + '.withorder')
	from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
	reader = TsvReader(infile, cnames = False)
	writer = TsvWriter(tmpfile)
	for r in reader:
		r[0] = mappings.get(r[0], r[0])
		writer.write(r)
	reader.close()
	writer.close()
	infile = tmpfile
	return dict(zip(mappings.values(), mappings.keys()))

def detach_order(mappings):
	if tool != 'sort' or not mappings:
		return
	shell.rm(infile)
	tmpfile = path.join(joboutdir, path.basename(infile) + '.withoutorder')
	from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
	reader = TsvReader(outfile, cnames = False)
	writer = TsvWriter(tmpfile)
	for r in reader:
		r[0] = mappings.get(r[0], r[0])
		writer.write(r)
	reader.close()
	writer.close()
	shell.mv(tmpfile, outfile)

tools = dict(
	sort     = run_sort,
	bedops   = run_bedops,
	bedtools = run_bedtools
)

try:
	if chrorder:
		mappings = attach_order()
	tools[tool]()
	if chrorder:
		detach_order(mappings)
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))
except:
	raise
finally:
	pass