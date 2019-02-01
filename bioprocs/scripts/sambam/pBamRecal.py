from sys import stderr
from pyppl import Box
from os import path
from bioprocs.utils import shell, mem2
from bioprocs.utils.poll import Poll
from bioprocs.utils.shell import Shell
from bioprocs.utils.reference import bamIndex

infile     = {{ i.infile | quote}}
bamIndex(infile, samtools = None)
outfile    = {{ o.outfile | quote}}
outprefix  = {{ o.outfile | prefix | quote}}
tool       = {{ args.tool | quote}}
gatk       = {{ args.gatk | quote}}
samtools   = {{ args.samtools | quote}}
picard     = {{ args.picard | quote}}
bamutil    = {{ args.bamutil | quote}}
params     = {{ args.params | repr}}
ref        = {{ args.ref | quote}}
tmpdir     = {{ args.tmpdir | quote}}
knownSites = {{ args.knownSites | quote}}
argsmem    = {{ args.mem | quote}}
nthread    = {{ args.nthread | quote}}
joboutdir  = {{ job.outdir | quote}}
joblen     = {{ proc.size | repr}}
jobindex   = {{ job.index | repr}}
workdir    = {{ proc.workdir | quote}}
tmpdir     = path.join (tmpdir, "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
# pass the index file
shell.ln_s(infile + '.bai', outfile + '.bai')
if not path.exists (tmpdir):
	shell.mkdir (tmpdir, p = True)
shell.TOOLS.update(dict(
	gatk     =  gatk,
	samtools = samtools,
	bamutil  = bamutil
))

# checks
if tool == 'gatk' and not path.isfile(knownSites):
	raise ValueError("knownSites file is required by GATK but is not specified (args.knownSites) or not exists.")

def run_gatk():
	mem = mem2(argsmem, '-jdict')
	intfile = path.join(joboutdir, outprefix + '.intervals')

	mem['-Djava.io.tmpdir={}'.format(shell.shquote(tmpdir))] = True
	gatksh = Shell(equal = ' ', dash = '-').gatk

	rtcparams    = params.get('RealignerTargetCreator', Box())
	rtcparams.T  = 'RealignerTargetCreator'
	rtcparams.R  = ref
	rtcparams.I  = infile
	rtcparams.o  = intfile
	rtcparams.nt = nthread
	rtcparams._  = list(mem.keys())
	gatksh(**rtcparams).run()

	bamfileir    = path.join(joboutdir, outprefix + '.ir.bam')
	irparams     = params.get('IndelRealigner', Box())
	irparams.T   = 'IndelRealigner'
	irparams.R   = ref
	irparams.I   = infile
	irparams.o   = bamfileir
	irparams._   = list(mem.keys())

	irparams.targetIntervals = intfile
	gatksh(**irparams).run()

	recaltable   = path.join(joboutdir, outprefix + '.recaltable')
	brparams     = params.get('BaseRecalibrator', Box())
	brparams.T   = 'BaseRecalibrator'
	brparams.R   = ref
	brparams.I   = bamfileir
	brparams.o   = recaltable
	brparams.nct = nthread
	brparams._   = list(mem.keys())

	brparams.knownSites = knownSites
	gatksh(**brparams).run()

	prparams     = params.get('PrintReads', Box())
	prparams.T   = 'PrintReads'
	prparams.R   = ref
	prparams.I   = bamfileir
	prparams.o   = outfile
	prparams.nct = nthread
	prparams._   = list(mem.keys())
	gatksh(**prparams).run()

	shell.rm(bamfileir, f = True)
	shell.mv(outprefix + '.bai', outfile + '.bai')

def run_bamutil():
	bush = Shell(subcmd = True, equal = ' ').bamutil
	if knownSites:
		params.dbsnp = knownSites
	params['in'] = infile

	params.out     = outfile
	params.refFile = ref
	params.verbose = True
	refcache       = path.splitext(ref)[0] + '-bs.umfa'
	if not path.isfile(refcache):
		poll = Poll(workdir, joblen, jobindex)
		poll.first(lambda **kwargs: bush.recal(**kwargs).run(), **params)
	else:
		bush.recab(**params).run()

tools = dict(
	gatk    = run_gatk,
	bamutil = run_bamutil
)
try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	shell.rm(tmpdir, r = True, f = True)
