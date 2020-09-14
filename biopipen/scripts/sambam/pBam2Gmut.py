from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell, mem2
from bioprocs.utils.reference import bamIndex

ref       = {{args.ref | quote}}
jobindex  = {{job.index | repr}}
tmpdir    = {{args.tmpdir | quote}}
procid    = {{proc.id | quote}}
infile    = {{i.infile | quote}}
prefix    = {{i.infile | fn2 | quote}}
tool      = {{args.tool | quote}}
gz        = {{args.gz | bool}}
outfile   = {{o.outfile | quote}}
params    = {{args.params | repr}}
samtools  = {{args.samtools | quote}}
mem       = {{args.mem | quote}}
nthread   = {{args.nthread | repr}}
gatk      = {{args.gatk | quote}}
cfgParams = {{args.cfgParams | repr}}
platypus  = {{args.platypus | quote}}
ssniffer  = {{args.snvsniffer | quote}}
strelka   = {{args.strelka | quote}}
vardict   = {{ args.vardict | quote}}
joboutdir = {{job.outdir | quote}}
bamIndex(infile,  samtools = None, nthread = nthread)
shell.load_config(dict(
	gatk     = gatk,
	platypus = platypus,
	ssniffer = ssniffer,
	strelka  = strelka,
	vardict  = vardict
))

if not path.exists(ref):
	raise ValueError("Reference file not specified")

tmpdir = path.join(tmpdir, "{procid}.{prefix}.{jobindex}".format(
	jobindex = jobindex,
	procid   = procid,
	prefix  = prefix
))
if not path.exists (tmpdir):
	shell.mkdir(tmpdir, p = True)

if gz:
	outfile = outfile[:-3]

def run_gatk():
	gatkmem = mem2(mem, 'jdict')
	gatkmem['Djava.io.tmpdir={!r}'.format(tmpdir)] = True

	params.T   = 'HaplotypeCaller'
	params.R   = ref
	params.I   = infile
	params.o   = outfile
	params.threads = nthread
	params._   = list(gatkmem.keys())

	shell.gatk(**params).fg
	if gz: shell.gzip(outfile)

def run_vardict():
	params.v      = True
	params.G      = ref
	params.b      = infile
	# params._out   = outfile
	# params._debug = True
	shell.vardict(**params).r > outfile
	if gz: shell.gzip(outfile)

def run_snvsniffer():
	hfile = path.join(joboutdir, prefix + '.header')
	shell.samtools.view(H = True).r > hfile

	params.g = ref
	params.o = outfile
	params._ = [hfile, infile]
	shell.ssniffer.snp(**params).fg
	if gz: shell.gzip(outfile)

def run_platypus():
	params.refFile     = ref
	params.bamFiles    = infile
	params.nCPU        = nthread
	params.output      = outfile
	params.logFileName = outfile + '.platypus.log'
	shell.platypus.callVariants(**params).fg
	if gz: shell.gzip(outfile)

def run_strelka():
	# config
	cfgParams.bam            = infile
	cfgParams.referenceFasta = ref
	cfgParams.runDir         = joboutdir
	shell.strelka(**cfgParams).fg

	# run the pipeline
	params.m = 'local'
	params.j = nthread
	params.g = int(float(mem2(mem, 'G')[:-1]))
	params._exe = path.join(joboutdir, 'runWorkflow.py')
	shell.runWorkflow(**params).fg

	# mv output file to desired outfile
	ofile = path.join(joboutdir, 'results', 'variants', 'genome.S1.vcf.gz')
	shell.mv(ofile, outfile + '.gz')
	if not gz: shell.gunzip(outfile + '.gz')


tools = {
	'gatk'      : run_gatk,
	'snvsniffer': run_snvsniffer,
	'strelka'   : run_strelka,
	'platypus'  : run_platypus,
	'vardict'   : run_vardict
}

try:
	tools[tool]()
except Exception:
	raise
finally:
	shell.rm(tmpdir, r = True, f = True)
