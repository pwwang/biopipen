from os import path, makedirs, rename
from shutil import rmtree
from sys import stderr, exit
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, mem2
from bioprocs.utils.reference import bamIndex

{% python from os import path %}

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
cfgParams = {{args.cfgParams | quote}}
platypus  = {{args.platypus | quote}}
ssniffer  = {{args.snvsniffer | quote}}
strelka   = {{args.strelka | quote}}
vardict   = {{ args.vardict | quote}}
joboutdir = {{job.outdir | quote}}
bamIndex(infile,  samtools = None, nthread = nthread)

if not path.exists (ref):
	stderr.write ("Reference file not specified")
	exit (1)

tmpdir = path.join(tmpdir, "{procid}.{prefix}.{jobindex}".format(
	jobindex = jobindex,
	procid   = procid,
	prefix  = prefix
))
if not path.exists (tmpdir):
	makedirs (tmpdir)

if gz:	
	outfile = outfile[:-3]

def run_gatk():
	params.R = ref
	params.I = infile
	params.o = outfile
	params.nct = nthread
	cmd = '{gatk} -T HaplotypeCaller {mem} -Djava.io.tmpdir={tmpdir!r} {args}'.format(
		gatk   = gatk,
		mem    = mem2(mem, 'java'),
		tmpdir = tmpdir,
		args   = cmdargs(params, dash = '-', equal = ' ')
	)
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_vardict():
	params.v = True
	params.G = ref
	params.b = infile
	cmd = '{vardict} {args} > {outfile!r}'.format(
		vardict = vardict,
		args    = cmdargs(params),
		outfile = outfile
	)
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_snvsniffer():
	hfile = path.join(joboutdir, prefix + '.header')
	runcmd('{samtools} view -H {infile!r} > {hfile!r}'.format(
		samtools = samtools,
		infile   = infile,
		hfile    = hfile
	))

	params.g   = ref
	params.o   = outfile
	params[''] = [hfile, infile]
	cmd = '{ssniffer} snp {args}'.format(ssniffer = ssniffer, args = cmdargs(params))
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_platypus():
	params.refFile     = ref
	params.bamFiles    = infile
	params.nCPU        = nthread
	params.output      = outfile
	params.logFileName = outfile + '.platypus.log'
	cmd = '{platypus} callVariants {args}'.format(platypus = platypus, args = cmdargs(params))
	runcmd(cmd)
	if gz: runcmd(['gzip', outfile])

def run_strelka():
	# config
	cfgParams.bam            = infile
	cfgParams.referenceFasta = ref
	cfgParams.runDir         = joboutdir
	runcmd('{strelka} {args}'.format(strelka = strelka, args = cmdargs(cfgParams)))

	# run the pipeline
	params.m = 'local'
	params.j = nthread
	params.g = mem2(mem, 'G')[:-1]
	cmd = '{runWorkflow} {args}'.format(
		runWorkflow = path.join(joboutdir, 'runWorkflow.py'),
		args        = cmdargs(params)
	)
	runcmd(cmd)

	# mv output file to desired outfile
	ofile = path.join(joboutdir, 'results', 'variants', 'genome.S1.vcf.gz')
	rename(ofile, outfile + '.gz')
	if not gz:
		runcmd(['gunzip', outfile + '.gz'])


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
	rmtree (tmpdir)

