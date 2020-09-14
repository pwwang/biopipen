"""Bam file recalibration"""
# pylint: disable=not-a-mapping,unsupported-assignment-operation
from sys import stderr
from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell, mem2
from bioprocs.utils.poll import Poll
from bioprocs.utils.reference import bam_index

infile     = {{ i.infile | quote}}
bam_index(infile, samtools=None)
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
# not work, have to re-index the output file
#shell.ln_s(infile + '.bai', outfile + '.bai')
if not path.exists (tmpdir):
    shell.mkdir (tmpdir, p = True)
shell.load_config(
    gatk     =  gatk,
    samtools = samtools,
    bamutil  = bamutil
)

# checks
if tool == 'gatk' and not path.isfile(knownSites):
    raise ValueError("knownSites file is required by GATK but "
                     "is not specified (args.knownSites) or not exists.")

def run_gatk():
    """Run GATK"""
    mem = mem2(argsmem, '-jdict')
    intfile = path.join(joboutdir, outprefix + '.intervals')

    mem['-Djava.io.tmpdir={}'.format(tmpdir)] = True

    rtcparams    = params.get('RealignerTargetCreator', Diot())
    rtcparams.T  = 'RealignerTargetCreator'
    rtcparams.R  = ref
    rtcparams.I  = infile
    rtcparams.o  = intfile
    rtcparams.nt = nthread
    rtcparams._  = list(mem.keys())
    shell.gatk(**rtcparams).fg

    bamfileir    = path.join(joboutdir, outprefix + '.ir.bam')
    irparams     = params.get('IndelRealigner', Diot())
    irparams.T   = 'IndelRealigner'
    irparams.R   = ref
    irparams.I   = infile
    irparams.o   = bamfileir
    irparams._   = list(mem.keys())

    irparams.targetIntervals = intfile
    shell.gatk(**irparams).fg

    recaltable   = path.join(joboutdir, outprefix + '.recaltable')
    brparams     = params.get('BaseRecalibrator', Diot())
    brparams.T   = 'BaseRecalibrator'
    brparams.R   = ref
    brparams.I   = bamfileir
    brparams.o   = recaltable
    brparams.nct = nthread
    brparams._   = list(mem.keys())

    brparams.knownSites = knownSites
    shell.gatk(**brparams).fg

    prparams     = params.get('PrintReads', Diot())
    prparams.T   = 'PrintReads'
    prparams.R   = ref
    prparams.I   = bamfileir
    prparams.o   = outfile
    prparams.nct = nthread
    prparams._   = list(mem.keys())
    shell.gatk(**prparams).fg

    shell.rm_rf(bamfileir)
    shell.mv(outprefix + '.bai', outfile + '.bai')

def run_bamutil():
    """Run bamutil"""
    if knownSites:
        params.dbsnp = knownSites
    params['in'] = infile

    params.out     = outfile
    params.refFile = ref
    params.verbose = True
    refcache       = path.splitext(ref)[0] + '-bs.umfa'
    if not path.isfile(refcache):
        poll = Poll(workdir, joblen, jobindex)
        poll.first(shell.bamutil.recab, **params)
    else:
        shell.bamutil.recab(**params).fg
    bam_index(outfile, samtools=samtools)

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
