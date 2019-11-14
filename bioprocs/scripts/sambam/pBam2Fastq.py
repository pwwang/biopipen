from os import makedirs, path
from pyppl import Box
from bioprocs.utils import logger, mem2, shell2 as shell

# bam2fastq will create {i.infile}.tmp, use file in indir in case of permission issue
infile  = {{i.infile | quote}}
fqfile1 = {{o.fqfile1 | quote}}
fqfile2 = {{o.fqfile2 | quote}}
tool    = {{args.tool | quote}}
params  = {{args.params | repr}}
mem     = {{ args.mem | quote }}
gz      = {{args.gz | bool}}
tmpdir  = {{args.tmpdir | quote}}

biobambam = {{args.biobambam | quote}}
bedtools  = {{args.bedtools | quote}}
samtools  = {{args.samtools | quote}}
picard    = {{args.picard | quote}}

shell.load_config(
	biobambam = biobambam,
	bedtools  = bedtools,
	samtools  = samtools,
	picard    = picard
)

if gz:
	fqfile1 = fqfile1[:-3]
	fqfile2 = fqfile2[:-3]

tmpdir = path.join(tmpdir, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

params._out_ = {{job.outfile | quote}}
params._err_ = {{job.errfile | quote}}

def gzip_fqs():
	if gz:
		shell.gzip(fqfile1)
		shell.gzip(fqfile2)

def run_biobambam():
	params.gz       = 0
	params.F        = fqfile1
	params.F2       = fqfile2
	params.T        = path.join(tmpdir, infile + '.tmp')
	params.filename = infile
	if infile.endswith('.sam'):
		params.inputformat = 'sam'
	shell.biobambam(**params)
	gzip_fqs()

def run_bedtools():
	params.i = infile
	params.fq = fqfile1
	params.fq2 = fqfile2
	shell.bedtools.bamtofastq(**params)
	gzip_fqs()

def run_samtools():
	params.t = True
	params['1'] = fqfile1
	params['2'] = fqfile2
	params._ = infile
	shell.samtools.fastq(**params)
	gzip_fqs()

def run_picard():
	mem = mem2({{ args.mem | quote }}, 'Java')
	params[mem]                = True
	params['-Djava.io.tmpdir'] = tmpdir
	params.TMP_DIR = tmpdir
	params.I = infile
	params.F = fqfile1
	params.F2 = fqfile2
	shell.picard(**params)
	gzip_fqs()

tools = dict(
	biobambam = run_biobambam,
	bedtools = run_bedtools,
	samtools = run_samtools,
	picard = run_picard,
)

try:
	tools[tool]()
except KeyError:
	raise ValueError('Tool "%r" not supported.' % tool)
finally:
	shell.rm_rf(tmpdir)
