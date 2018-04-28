
from os import path, remove, makedirs
from sys import stderr, exit
from time import sleep
from datetime import date

from pyppl import Box
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.poll import Poll

gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{out.outfile | quote}}
if gz:	outfile = outfile[:-3]

poll = Poll({{proc.workdir | quote}}, {{proc.size}}, {{job.index}})

params = {{args.params}}

############## cnvkit
{% if args.tool | lambda x: x=='cnvkit' %}

targetDone    = "{{proc.workdir}}/1/output/cnvkit_target.done"
referenceDone = "{{proc.workdir}}/1/output/cnvkit_reference.done"
targetCov     = "{{out.outdir}}/{{in.infile | fn}}.targetcov.cnn"
accessfile    = '{{proc.workdir}}/1/output/cnvkit_access.bed'
targetfile    = '{{proc.workdir}}/1/output/cnvkit_targets.bed'
refcnn        = '{{proc.workdir}}/1/output/reference.cnn'
fixedCnr      = "{{out.outdir}}/{{in.infile | fn}}.cnr"
segfile       = "{{out.outdir}}/{{in.infile | fn}}.cns"
callfile      = "{{out.outdir}}/{{in.infile | fn}}.call.cns"
# report files
breaksfile    = "{{out.outdir}}/{{in.infile | fn}}.breaks.txt"
gainlossfile  = "{{out.outdir}}/{{in.infile | fn}}.gainloss.txt"
metricsfile   = "{{out.outdir}}/{{in.infile | fn}}.metrics.txt"
segmetricsfile= "{{out.outdir}}/{{in.infile | fn}}.segmetrics.txt"
openblas_nthr = "export OPENBLAS_NUM_THREADS={{args.nthread}}; export OMP_NUM_THREADS={{args.nthread}}; export NUMEXPR_NUM_THREADS={{args.nthread}}; export MKL_NUM_THREADS={{args.nthread}}; "

cnvkitAccessParams = {{args.cnvkitAccessParams}}
cnvkitAccessParams['o'] = accessfile
cmd1 = '%s {{args.cnvkit}} access "{{args.ref}}" %s' % (openblas_nthr, cmdargs(cnvkitAccessParams))

cnvkitTargetParams = {{args.cnvkitTargetParams}}
cnvkitTargetParams['o'] = targetfile
cmd2 = '{{args.cnvkit}} target "%s" %s' % (accessfile, cmdargs(cnvkitTargetParams))

poll.first(cmd1 + '; ' + cmd2)

#pollingNon1st ({{job.index}}, cmd1 + '; ' + cmd2, targetDone, t = 60)

cnvkitCoverageParams = {{args.cnvkitCoverageParams}}
cnvkitCoverageParams['p'] = {{args.nthread}}
cnvkitCoverageParams['o'] = targetCov
cmd = '%s {{args.cnvkit}} coverage "{{in.infile}}" "%s" %s' % (openblas_nthr, targetfile, cmdargs(cnvkitCoverageParams))

poll.all(cmd)
#pollingAll ({{proc.workdir | quote}}, {{proc.size}}, {{job.index}}, cmd, "cnvkit_coverage.done")

cnvkitReferenceParams = {{args.cnvkitReferenceParams}}
cnvkitReferenceParams['o'] = refcnn
cnvkitReferenceParams['f'] = {{args.ref | quote}}
cmd = '%s {{args.cnvkit}} reference {{proc.workdir}}/*/output/*/*.targetcov.cnn %s' % (openblas_nthr, cmdargs(cnvkitReferenceParams))

poll.first(cmd)
#pollingNon1st ({{job.index}}, cmd, referenceDone)

mtfile = "{{out.outdir}}/cnvkit_mt"
open(mtfile, 'w').close()

cnvkitFixParams = {{args.cnvkitFixParams}}
cnvkitFixParams['o'] = fixedCnr
cmd = '%s {{args.cnvkit}} fix "%s" "%s" "%s" %s' % (openblas_nthr, targetCov, mtfile, refcnn, cmdargs(cnvkitFixParams))
runcmd (cmd)

if path.getsize(fixedCnr) < 60:
	open(segfile, 'w').write('chromosome	start	end	gene	log2	depth	probes	weight\\n')
else:
	cnvkitSegmentParams = {{args.cnvkitSegmentParams}}
	cnvkitSegmentParams['o'] = segfile
	cnvkitSegmentParams['p'] = {{args.nthread}}

	cmd = '%s {{args.cnvkit}} segment %s "%s"' % (openblas_nthr, cmdargs(cnvkitSegmentParams), fixedCnr)
	runcmd (cmd)

if path.getsize(segfile) < 60:
	open(callfile, 'w').write('chromosome	start	end	gene	log2	cn	depth	probes	weight\\n')
else:
	cnvkitCallParams = {{args.cnvkitCallParams}}
	cnvkitCallParams['o'] = callfile
	cmd = '%s {{args.cnvkit}} call %s "%s"' % (openblas_nthr, cmdargs(cnvkitCallParams), segfile)
	runcmd (cmd)

# reports
{% if args.cnvkitReport %}
cnvkitBreaksParams = {{args.cnvkitBreaksParams}}
cnvkitBreaksParams['o'] = breaksfile
cmd = '%s {{args.cnvkit}} breaks "%s" "%s" %s' % (openblas_nthr, fixedCnr, callfile, cmdargs(cnvkitBreaksParams))
runcmd (cmd)

cnvkitGainlossParams = {{args.cnvkitGainlossParams}}
cnvkitGainlossParams['s'] = callfile
cnvkitGainlossParams['o'] = gainlossfile
cmd = '%s {{args.cnvkit}} gainloss "%s" %s' % (openblas_nthr, fixedCnr, cmdargs(cnvkitGainlossParams))
runcmd (cmd)

cnvkitMetricsParams = {{args.cnvkitMetricsParams}}
cnvkitMetricsParams['s'] = callfile
cnvkitMetricsParams['o'] = metricsfile
cmd = '%s {{args.cnvkit}} metrics "%s" %s' % (openblas_nthr, fixedCnr, cmdargs(cnvkitMetricsParams))
runcmd (cmd)

cnvkitSegmetricsParams = {{args.cnvkitSegmetricsParams}}
cnvkitSegmetricsParams['s'] = callfile
cnvkitSegmetricsParams['o'] = segmetricsfile
cmd = '%s {{args.cnvkit}} segmetrics "%s" %s' % (openblas_nthr, fixedCnr, cmdargs(cnvkitSegmetricsParams))
runcmd (cmd)
{% endif %}

# plots
{% if args.cnvkitPlot %}
{% for i, param in args.cnvkitScatterParams | lambda x: enumerate(x) %}
param = {{param}}
param['s'] = callfile
param['o'] = "{{out.outdir}}/{{in.infile | fn}}.scatter{{i | lambda x: x+1}}.pdf"
cmd = '%s {{args.cnvkit}} scatter "%s" %s' % (openblas_nthr, fixedCnr, cmdargs(param))
runcmd (cmd)
{% endfor %}

cnvkitDiagramParams = {{args.cnvkitDiagramParams}}
cnvkitDiagramParams['s'] = callfile
cnvkitDiagramParams['o'] = "{{out.outdir}}/{{in.infile | fn}}.diagram.pdf"
cmd = '%s {{args.cnvkit}} diagram "%s" %s' % (openblas_nthr, fixedCnr, cmdargs(cnvkitDiagramParams))
runcmd (cmd)

poll.all(lambda x: not path.exists(x), callfile)
#pollingAll ({{proc.workdir | quote}}, {{proc.size}}, {{job.index}}, 'if [[ ! -e "%s" ]]; then sleep 1; fi' % callfile, "cnvkit_call.done")

{% if job.index | lambda x: x == 0 %}
{% for i, param in args.cnvkitHeatmapParams | lambda x: enumerate(x) %}
param = {{param}}
param['o'] = "{{out.outdir}}/{{in.infile | fn}}.heatmap{{i | lambda x: x+1}}.pdf"
cmd = '%s {{args.cnvkit}} heatmap "{{proc.workdir}}"/*/output/*/*.call.cns %s' % (openblas_nthr, cmdargs(param))
runcmd (cmd)
{% endfor %}
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
		fout.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{{in.infile | fn}}\\n')
else:
	cnvkitExportParams = {{args.cnvkitExportParams}}
	cnvkitExportParams['o'] = outfile
	cmd = '%s {{args.cnvkit}} export vcf "%s" %s' % (openblas_nthr, callfile, cmdargs(cnvkitExportParams))
	runcmd (cmd)

############## cnvnator
{% elif args.tool | lambda x: x=='cnvnator' %}
rootfile = '{{out.outdir}}/{{in.infile | fn}}.root'
callfile = '{{out.outdir}}/{{in.infile | fn}}.cnvnator'
cmd = '{{args.cnvnator}} -root "%s" -genome {{args.cnvnatorGenome}} -tree "{{in.infile}}" ' % (rootfile)
runcmd (cmd)
cmd = '{{args.cnvnator}} -root "%s" -genome {{args.cnvnatorGenome}} -his {{args.cnvnatorBinsize}}' % (rootfile)
runcmd (cmd)
cmd = '{{args.cnvnator}} -root "%s" -stat {{args.cnvnatorBinsize}}' % (rootfile)
runcmd (cmd)
cmd = '{{args.cnvnator}} -root "%s" -partition {{args.cnvnatorBinsize}}' % (rootfile)
runcmd (cmd)
cmd = '{{args.cnvnator}} -root "%s" -call {{args.cnvnatorBinsize}} > "%s"' % (rootfile, callfile)
runcmd (cmd)
cmd = 'cd "{{out.outdir}}"; {{args.cnvnator2vcf}} "%s" > "%s"' % (path.basename(callfile), outfile)
runcmd (cmd)

############## wandy
{% elif args.tool | lambda x: x=='wandy' %}
# get Wandy tool.info
from distutils.spawn import find_executable
toolinfo   = path.join (path.dirname(find_executable("{{args.wandy}}")), "tool.info")
retdir     = "{{out.outdir}}"
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

params['i'] = {{in.infile | quote}}
params['t'] = {{args.wandyType}}
cmd = 'cd "%s"; {{args.wandy}} %s' % (retdir, cmdargs(params))
runcmd (cmd)

# TODO: convert to vcf
open({{out.outfile | quote}}, 'w').close()
{% endif %}

if gz:	runcmd ('gzip "%s"' % (outfile))
