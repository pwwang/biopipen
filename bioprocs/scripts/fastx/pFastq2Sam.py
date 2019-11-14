import re
from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

# variables
infile1   = {{i.fq1 | quote}}
infile2   = {{i.fq2 | quote}}
outfile   = {{o.outfile | quote}}
outprefix = {{o.outfile | fn | quote}}
outfmt    = {{args.outfmt | quote}}
outdir    = {{job.outdir | quote}}
ref       = {{args.ref | quote}}
rg        = {{args.rg | repr}}
params    = {{args.params | repr}}
samtools  = {{args.samtools | quote}}
bowtie2   = {{args.bowtie2 | quote}}
bwa       = {{args.bwa | quote}}
ngm       = {{args.ngm | quote}}
star      = {{args.star | quote}}
tool      = {{args.tool | quote}}
nthread   = {{args.nthread}}
jobindex  = {{job.index}}

shell.load_config(dict(
	bowtie2  = bowtie2,
	samtools = samtools,
	bwa      = bwa,
	star     = star,
	ngm      = ngm,
))

# detemine default read group
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', outprefix)
	rg['ID'] = g.group(1) if g else "{outprefix}.L{jobindex}".format(outprefix = outprefix, jobindex = jobindex)
if not rg['SM']:
	rg['SM'] = outprefix

def sam2bam(samfile, bamfile):
	shell.samtools.view(
		S = True, b = True, _out = bamfile, _ = samfile, **{'@': nthread}
	)
	shell.rm_rf(samfile)

def run_bowtie2():
	params.threads = nthread
	params.x = ref
	params.S = outfile if outfmt == 'sam' else path.splitext(outfile)[0] + '.sam'
	params['1'] = infile1
	params['2'] = infile2
	params['rg-id'] = rg['ID']
	params['rg'] = ['{}:{}'.format(k, v) for k, v in rg.items() if k != 'ID']
	shell.fg.bowtie2(**params)
	if outfmt == 'bam': sam2bam(params.S, outfile)

def run_bwa():
	params.t = nthread
	params.R = "@RG\\tID:{id}\\t{rg}".format(id = rg['ID'], rg = "\\t".join(
		'{}:{}'.format(k, v) for k, v in rg.items() if k != 'ID'
	))
	params.o = outfile if outfmt == 'sam' else path.splitext(outfile)[0] + '.sam'
	params._ = [ref, infile1, infile2]
	shell.fg.bwa.mem(**params)
	if outfmt == 'bam': sam2bam(params.o, outfile)

def run_ngm():
	params['1'] = infile1
	params['2'] = infile2
	params.r = ref
	params.o = outfile
	params.b = outfmt == 'bam'
	params.t = nthread
	for k, v in rg.items():
		params['rg-' + k.lower()] = v
	shell.fg.ngm(**params)

def run_star():
	params.genomeDir        = path.splitext(ref)[0] + '.star'
	params.readFilesIn      = [infile1, infile2]
	params.readFilesCommand = ("cat", "zcat", "bzcat")[
		1 if infile1.endswith('.gz') else 2 if infile1.endswith('.bz2') else 0
	]
	params.readNameSeparator = '.'
	params.outFileNamePrefix = outdir + '/'
	params.outSAMtype        = [outfmt.upper(), 'Unsorted']
	shell.fg.star(**params)

	starout = path.join(outdir, "Aligned.out.{}".format(outfmt))
	if path.isfile(starout):
		shell.mv(starout, outfile)

tools = dict(
	bowtie2 = run_bowtie2,
	bwa     = run_bwa,
	ngm     = run_ngm,
	star    = run_star
)

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))

