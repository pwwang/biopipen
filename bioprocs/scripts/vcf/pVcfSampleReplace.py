from os import path
from pyppl.utils import alwaysList
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
samfile  = {{i.samfile | quote}}
outfile  = {{o.outfile | quote}}
samples  = {{args.samples | ?:isinstance(_, str) and _.startswith('lambda') | !repr}}
nthread  = {{args.nthread | quote}}
bcftools = {{args.bcftools | quote}}

if isinstance(samples, str):
	if path.isfile(samples):
		with open(samples) as f:
			samples = [line.strip() for line in f.readlines() if line.strip()]
	else:
		samples = alwaysList(samples)
elif samples and not callable(samples):
	samples = list(samples)

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
shell.bcftools.reheader(**params).fg
