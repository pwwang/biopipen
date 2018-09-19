
from os import path, remove, makedirs
from sys import stderr, exit
from time import sleep
from datetime import date

from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, log2pyppl, logger
from bioprocs.utils.poll import Poll

gz       = {{args.gz | lambda x: bool(x)}}
infile   = {{i.infile | quote}}
infn     = {{i.infile | fn | quote}}
outfile  = {{o.outfile | quote}}
outdir   = {{o.outdir | quote}}
workdir  = {{proc.workdir | quote}}
nthread  = {{args.nthread | repr}}
cnvnator = {{args.cnvnator | repr}}
cnvkit   = {{args.cnvkit | quote}}
wandy    = {{args.wandy | quote}}
ref      = {{args.ref | quote}}
genome   = {{args.genome | quote}}
params   = {{args.params}}

def log2log(*args, **kwargs):
	# debug = False
	# pass 
	logger.info(*args, **kwargs)
	# log2pyppl(*args, **kwargs)

cnvkitParamsKeys = [
	'access', 'target', 'coverage', 'reference', 'fix', 'segment', 'call', 'plot', 'breaks',
	'gainloss', 'metrics', 'segmetrics', 'export', 'scatter', 'heatmap', 'diagram'
]

if gz:	outfile = outfile[:-3]

poll = Poll(workdir, {{proc.size}}, {{job.index}})

############## cnvkit
{% if args.tool == 'cnvkit' %}

#log2log('CNVkit: Initialize ...')

targetCov     = "{outdir}/{infn}.targetcov.cnn".format(outdir = outdir, infn = infn)
accessfile    = "{workdir}/1/output/cnvkit_access.bed".format(workdir = workdir)
targetfile    = "{workdir}/1/output/cnvkit_targets.bed".format(workdir = workdir)
refcnn        = "{workdir}/1/output/reference.cnn".format(workdir = workdir)
fixedCnr      = "{outdir}/{infn}.cnr".format(outdir = outdir, infn = infn)
segfile       = "{outdir}/{infn}.cns".format(outdir = outdir, infn = infn)
callfile      = "{outdir}/{infn}.call.cns".format(outdir = outdir, infn = infn)
# report files
breaksfile    = "{outdir}/{infn}.breaks.txt".format(outdir = outdir, infn = infn)
gainlossfile  = "{outdir}/{infn}.gainloss.txt".format(outdir = outdir, infn = infn)
metricsfile   = "{outdir}/{infn}.metrics.txt".format(outdir = outdir, infn = infn)
segmetricsfile= "{outdir}/{infn}.segmetrics.txt".format(outdir = outdir, infn = infn)
openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}; ".format(nthread = nthread)

cnvkitAccessParams      = params.access
cnvkitAccessParams['o'] = accessfile
cmd1 = openblas_nthr + '{cnvkit} access {ref} {args}'.format(cnvkit = cnvkit, ref = repr(ref), args = cmdargs(cnvkitAccessParams))

cnvkitTargetParams = params.target
cnvkitTargetParams['o'] = targetfile
cmd2 = '{cnvkit} target {accessfile} {args}'.format(cnvkit = cnvkit, accessfile = repr(accessfile), args = cmdargs(cnvkitTargetParams))

log2log('CNVkit: Run access and target at job #0 ...')
poll.first(cmd1 + '; ' + cmd2, lockfile = 'access.poll.lock')
log2log('CNVkit: Run access and target at job #0 ... done')

cnvkitCoverageParams = params.coverage
cnvkitCoverageParams['p'] = nthread
cnvkitCoverageParams['o'] = targetCov
cmd = openblas_nthr + '{cnvkit} coverage {infile} {targetfile} {args}'.format(cnvkit = cnvkit, infile = repr(infile), targetfile = repr(targetfile), args = cmdargs(cnvkitCoverageParams))

log2log('CNVkit: Run coverage at all jobs ...')
poll.all(cmd, lockfile = 'coverage.poll.lock')
log2log('CNVkit: Run coverage at all jobs ... done')

cnvkitReferenceParams = params.reference
cnvkitReferenceParams['o'] = refcnn
cnvkitReferenceParams['f'] = ref
cmd = openblas_nthr + '{cnvkit} reference {workdir}/*/output/*/*.targetcov.cnn {args}'.format(cnvkit = cnvkit, workdir = repr(workdir), args = cmdargs(cnvkitReferenceParams))

log2log('CNVkit: Run reference at job #0 ...')
poll.first(cmd, lockfile = 'reference.poll.lock')
log2log('CNVkit: Run reference at job #0 ... done')

mtfile = "{outdir}/cnvkit_mt".format(outdir = outdir)
open(mtfile, 'w').close()

cnvkitFixParams      = params.fix
cnvkitFixParams['o'] = fixedCnr
cmd = openblas_nthr + '{cnvkit} fix {targetCov} {mtfile} {refcnn} {args}'.format(cnvkit = cnvkit, targetCov = repr(targetCov), mtfile = repr(mtfile), refcnn = repr(refcnn), args = cmdargs(cnvkitFixParams))

log2log('CNVkit: Run fix at all jobs ...')
runcmd (cmd)
log2log('CNVkit: Run fix at all jobs ... done')

if path.getsize(fixedCnr) < 60:
	open(segfile, 'w').write('chromosome	start	end	gene	log2	depth	probes	weight\\n')
else:
	cnvkitSegmentParams = params.segment
	cnvkitSegmentParams['o'] = segfile
	cnvkitSegmentParams['p'] = nthread

	cmd = openblas_nthr + '{cnvkit} segment {args} {fixedCnr}'.format(cnvkit = cnvkit, args = cmdargs(cnvkitSegmentParams), fixedCnr = repr(fixedCnr))
	log2log('CNVkit: Run segment at all jobs ...')
	runcmd (cmd)
	log2log('CNVkit: Run segment at all jobs ... done')

if path.getsize(segfile) < 60:
	open(callfile, 'w').write('chromosome	start	end	gene	log2	cn	depth	probes	weight\\n')
else:
	cnvkitCallParams = params.call
	cnvkitCallParams['o'] = callfile
	cmd = openblas_nthr + '{cnvkit} call {args} {segfile}'.format(cnvkit = cnvkit, args = cmdargs(cnvkitCallParams), segfile = repr(segfile))
	log2log('CNVkit: Run call at all jobs ...')
	runcmd (cmd)
	log2log('CNVkit: Run call at all jobs ... done')

# reports
{% if args.report %}
cnvkitBreaksParams = params.breaks
cnvkitBreaksParams['o'] = breaksfile
cmd = openblas_nthr + '{cnvkit} breaks {fixedCnr} {callfile} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), callfile = repr(callfile), args = cmdargs(cnvkitBreaksParams))
runcmd (cmd)

cnvkitGainlossParams = params.gainloss
cnvkitGainlossParams['s'] = callfile
cnvkitGainlossParams['o'] = gainlossfile
cmd = openblas_nthr + '{cnvkit} gainloss {fixedCnr} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), args = cmdargs(cnvkitGainlossParams))
runcmd (cmd)

cnvkitMetricsParams = params.metrics
cnvkitMetricsParams['s'] = callfile
cnvkitMetricsParams['o'] = metricsfile
cmd = openblas_nthr + '{cnvkit} metrics {fixedCnr} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), args = cmdargs(cnvkitMetricsParams))
runcmd (cmd)

cnvkitSegmetricsParams = params.segmetrics
cnvkitSegmetricsParams['s'] = callfile
cnvkitSegmetricsParams['o'] = segmetricsfile
cmd = openblas_nthr + '{cnvkit} segmetrics {fixedCnr} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), args = cmdargs(cnvkitSegmetricsParams))
runcmd (cmd)
{% endif %}

# plots
{% if args.plot %}
# require job.index == 0?
if not isinstance(params.scatter, list):
	params.scatter = [params.scatter]
for i, param in enumerate(params.scatter):
	param['s'] = callfile
	param['o'] = "{outdir}/{infn}.scatter{i}.pdf".format(outdir = outdir, infn = infn, i = i+1)
	cmd = openblas_nthr + '{cnvkit} scatter {fixedCnr} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), args = cmdargs(param))
	runcmd (cmd)

cnvkitDiagramParams = params.diagram
cnvkitDiagramParams['s'] = callfile
cnvkitDiagramParams['o'] = "{outdir}/{infn}.diagram.pdf".format(outdir = outdir, infn = infn)
cmd = openblas_nthr + '{cnvkit} diagram {fixedCnr} {args}'.format(cnvkit = cnvkit, fixedCnr = repr(fixedCnr), args = cmdargs(cnvkitDiagramParams))
runcmd (cmd)

poll.all(lambda x: not path.exists(x), callfile, lockfile = 'callfile.poll.lock')

{% if job.index == 0 %}
if not isinstance(params.heatmap, list):
	params.heatmap = [params.heatmap]
for i, param in enumerate(params.heatmap):
	param['o'] = "{outdir}/{infn}.heatmap{i}.pdf".format(outdir = outdir, infn = infn, i = i+1)
	cmd = openblas_nthr + '{cnvkit} heatmap {workdir}/*/output/*/*.call.cns {args}'.format(cnvkit = cnvkit, workdir = repr(workdir), args = cmdargs(param))
	runcmd (cmd)
{% endif %}
{% endif %}

# to vcf
if path.getsize(callfile) < 60: # no data generated
	with open (outfile, 'w') as fout:
		fout.write('##fileformat=VCFv4.2\\n')
		fout.write('##fileDate=%s\\n' % str(date.today()).replace('-', ''))
		fout.write('##source=CNVkit\\n')
		fout.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\\n')
		fout.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\\n')
		fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\\n')
		fout.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\\n')
		fout.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\\n')
		fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n')
		fout.write('##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description="Fold change">\\n')
		fout.write('##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">\\n')
		fout.write('##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">\\n')
		fout.write('##ALT=<ID=DEL,Description="Deletion">\\n')
		fout.write('##ALT=<ID=DUP,Description="Duplication">\\n')
		fout.write('##ALT=<ID=CNV,Description="Copy number variable region">\\n')
		fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n')
		fout.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">\\n')
		fout.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\\n')
		fout.write('##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">\\n')
		fout.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{infn}\\n'.format(infn = infn))
else:
	cnvkitExportParams = params.export
	cnvkitExportParams['o'] = outfile
	cmd = openblas_nthr + '{cnvkit} export vcf {callfile} {args}'.format(cnvkit = cnvkit, callfile = repr(callfile), args = cmdargs(cnvkitExportParams))
	runcmd (cmd)

############## cnvnator
{% elif args.tool == 'cnvnator' %}
rootfile = '{outdir}/{infn}.root'.format(outdir = outdir, infn = infn)
callfile = '{outdir}/{infn}}.cnvnator'.format(outdir = outdir, infn = infn)
binsize  = {{args.binsize | repr}}

cmd = '{cnvnator} -root {rootfile} -genome {genome} -tree {infile}'.format(cnvnator = cnvnator, rootfile = repr(rootfile), genome = genome, infile = repr(infile))
runcmd (cmd)

cmd = '{cnvnator} -root {rootfile} -genome {genome} -his {binsize}'.format(cnvnator = cnvnator, rootfile = repr(rootfile), genome = genome, binsize = binsize)
runcmd (cmd)

cmd = '{cnvnator} -root {rootfile} -stat {binsize}'.format(cnvnator = cnvnator, rootfile = repr(rootfile), binsize = binsize)
runcmd (cmd)

cmd = '{cnvnator} -root {rootfile} -partition {binsize}'.format(cnvnator = cnvnator, rootfile = repr(rootfile), binsize = binsize)
runcmd (cmd)

cmd = '{cnvnator} -root {rootfile} -call {binsize}'.format(cnvnator = cnvnator, rootfile = repr(rootfile), binsize = binsize)
runcmd (cmd)

cmd = 'cd {outdir}; {cnvnator2vcf} {callfile} > {outfile}'.format(outdir = repr(outdir), cnvnator2vcf = cnvnator2vcf, callfile = repr(path.basename(callfile)), outfile = repr(outfile))
runcmd (cmd)

############## wandy
{% elif args.tool == 'wandy' %}
# get Wandy tool.info
from distutils.spawn import find_executable
toolinfo   = path.join (path.dirname(find_executable(wandy)), "tool.info")
retdir     = outdir
#if not path.exists(retdir):
#	makedirs(retdir)
myToolinfo = path.join (retdir, "tool.info")
if not path.exists(toolinfo):
	stderr.write('Cannot find tool.info in wandy source directory.')
	exit(1)
with open (toolinfo) as fin, open(myToolinfo, 'w') as fout:
	toolinfos  = [line.strip() for line in fin if not line.strip().startswith('WANDY_MEM') and not line.strip().startswith('WANDY_PARALLELED_MODE')]
	toolinfos.append ('WANDY_MEM={{args.mem}}')
	toolinfos.append ('WANDY_PARALLELED_MODE=0')
	fout.write ("\n".join(toolinfos) + "\n")

params      = {key: val for key, val in  params if key not in cnvkitParamsKeys}
params['i'] = infile
params['t'] = {{args.type | repr}}
cmd = 'cd {retdir}; {wandy} {args}'.format(retdir = repr(retdir), wandy = wandy, args = cmdargs(params))
runcmd (cmd)

# TODO: convert to vcf
open(outfile, 'w').close()
{% endif %}

if gz:	runcmd ('gzip "%s"' % (outfile))
