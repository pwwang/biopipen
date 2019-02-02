from pyppl import Box
from bioprocs.utils import shell

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote }}
tool     = {{ args.tool | quote}}
bedtools = {{ args.bedtools | quote}}
bedops   = {{ args.bedops | quote}}
argsmem  = {{ args.mem | quote}}
sortby   = {{ args.by | quote}}
unique   = {{ args.unique | bool}}
params   = {{ args.params | repr}}
tmpdir   = {{ args.tmpdir | quote}}

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
	shell.Shell().grep(v = '^#', infile).pipe(dash = '-', equal = '=').sort(**params, __stdout = outfile).run()

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

tools = dict(
	sort     = run_sort,
	bedops   = run_bedops,
	bedtools = run_bedtools
)

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))
except:
	raise
finally:
	pass