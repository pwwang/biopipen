from os import path
from diot import Diot
from pyppl.utils import always_list
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
samfile  = {{i.samfile | quote}}
outfile  = {{o.outfile | quote}}
samples  = {{args.samples | ?:isinstance(_, str) and _.startswith('lambda') | !repr}}
params   = {{args.params | repr}}
keep     = {{args.keep | quote}}
bcftools = {{args.bcftools | quote}}

if isinstance(samples, str):
	if path.isfile(samples):
		with open(samples) as f:
			samples = [line.strip() for line in f.readlines() if line.strip()]
	else:
		samples = always_list(samples)
elif samples and not callable(samples):
	samples = list(samples)

shell.load_config(bcftools = bcftools)

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
	osams = shell.bcftools.query(l = infile).splitlines()
	osams = [s.strip() for s in osams if s.strip()]
	if callable(samples):
		# get samples
		rsams = [s for s in osams if s and samples(s)]
	elif isinstance(samples[0], int):
		rsams = [osams[s] for s in samples]
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

shell.bcftools.view(**params).fg
