from os import path
from pyppl import Box
from bioprocs.utils import shell

{% python from os import path %}
{% python from bioprocs.utils import alwaysList %}
{% assign samcalls = dict(list = list, readlines = readlines, alwaysList = alwaysList, func = lambda x: x) %}
{% assign samtypes = lambda x, path = path: 'none' if x is None else 'list' if isinstance(x, (list, tuple)) else 'readlines' if path.isfile(x) else 'func' if x.startswith('lambda') else 'alwaysList' %}
infile   = {{i.infile | quote}}
samfile  = {{i.samfile | quote}}
outfile  = {{o.outfile | quote}}
samples  = {{args.samples | samcalls.get(samtypes(args.samples), repr) }}
params   = {{args.params | repr}}
keep     = {{args.keep | quote}}
bcftools = {{args.bcftools | quote}}

shell.TOOLS.bcftools = bcftools
bcftools = shell.Shell(subcmd = True, equal = ' ').bcftools

if not samfile and not samples:
	raise ValueError('Require either `i.samfile` or `args.samples`')

if samfile and samples and not callable(samples):
	raise ValueError('Both `i.samfile` and `args.samples` provided, I dont know which one to use.')

if samfile:
	with open(samfile) as f:
		rsams = [s.strip() for s in f if s.strip()]
	if callable(samples):
		rsams = [s for s in rsams if samples(s)]
else:
	if callable(samples):
		# get samples 
		osams = bcftools.query(l = infile).run(save = 'stdout').stdout.splitlines()
		osams = [s.strip() for s in osams if s.strip()]
		rsams = [s for s in osams if s.strip() and samples(s.strip())]
	else:
		rsams = list(samples)

samfile = outfile + '.samples'
with open(samfile, 'w') as f:
	f.write(''.join(s + '\n' for s in rsams))

params.update(dict(
	S = samfile if keep else '^' + samfile,
	o       = outfile,
	_       = infile
))

bcftools.view(**params).run()