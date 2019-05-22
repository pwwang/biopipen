from bioprocs.utils import shell2 as shell

{% python from os import path %}
{% python from bioprocs.utils import alwaysList %}
{% assign samcalls = dict(list = list, readlines = readlines, alwaysList = alwaysList, func = lambda x: x) %}
{% assign samtypes = lambda x, path = path: 'none' if x is None else 'list' if isinstance(x, (list, tuple)) else 'readlines' if path.isfile(x) else 'func' if x.startswith('lambda') else 'alwaysList' %}
infile   = {{i.infile | quote}}
samfile  = {{i.samfile | quote}}
outfile  = {{o.outfile | quote}}
samples  = {{args.samples | samcalls.get(samtypes(args.samples), repr) }}
nthread  = {{args.nthread | quote}}
bcftools = {{args.bcftools | quote}}

shell.load_config(bcftools = bcftools)

if not samfile and not samples:
	raise ValueError('Require either `i.samfile` or `args.samples`')

if samfile and samples and not callable(samples):
	raise ValueError('Both `i.samfile` and `args.samples` provided, I dont know which one to use.')

osams = shell.bcftools.query(l = infile).splitlines()
osams = [s.strip() for s in osams if s.strip()]
if samfile:
	with open(samfile) as f:
		rsams = [s.strip() for s in f if s.strip()]
	if callable(samples):
		rsams = [samples(s) for s in rsams]
else:
	if callable(samples):
		# get samples
		rsams = [samples(s.strip()) for s in osams if s.strip()]
	else:
		rsams = list(samples)

if len(osams) != len(rsams):
	raise ValueError('Unequal # samples: {} in VCF, {} to replace'.format(len(osams), len(rsams)))

samfile = outfile + '.newsamples'
with open(samfile, 'w') as f:
	f.write(''.join(s + '\n' for s in rsams))

params = dict(
	threads = nthread,
	o       = outfile,
	_       = infile,
	s       = samfile
)
shell.fg.bcftools.reheader(**params)
