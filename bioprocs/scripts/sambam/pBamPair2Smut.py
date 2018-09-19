from os import path, makedirs
from shutil import rmtree
from sys import stderr
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2

if not path.exists ({{args.ref | quote}}):
	stderr.write ("Reference file not specified")
	exit (1)

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{i.tumor | fn | fn}}.{{i.normal | fn | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

ref     = {{args.ref | quote}}
tool    = {{args.tool | quote}}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{o.outfile | quote}}
params  = {{args.params}}
if gz:	outfile = outfile[:-3]
try:
{% case args.tool %}
	############## gatk
	{% if when 'gatk' %}

	intvfile = "{{job.outdir}}/interval.list"
	cmd      = '{{args.samtools}} idxstats "{{i.tumor}}" | head -n -1 | cut -f1 > "%s"' % (intvfile)
	runcmd (cmd)

	mem = mem2 ({{args.mem | quote}}, 'java')

	params['R']        = ref
	params['I:tumor']  = {{i.tumor | quote}}
	params['I:normal'] = {{i.normal | quote}}
	params['o']        = outfile
	params['nct']      = {{args.nthread}}
	params['L']        = intvfile

	cmd = '{{args.gatk}} -T MuTect2 %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '-', equal = ' '))
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## somaticsniper
	{% when 'somaticsniper' %}
	params['f'] = ref
	params['F'] = 'vcf'

	cmd = '{{args.somaticsniper}} %s "{{i.tumor}}" "{{i.normal}}" "%s"' % (cmdargs(params), outfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## snvsniffer
	{% when 'snvsniffer' %}
	theader = "{{job.outdir}}/{{i.tumor | bn}}.header"
	cmd = '{{args.samtools}} view -H "{{i.tumor}}" > "%s"' % theader
	runcmd (cmd)
	nheader = "{{job.outdir}}/{{i.normal | bn}}.header"
	cmd = '{{args.samtools}} view -H "{{i.normal}}" > "%s"' % nheader
	runcmd (cmd)

	params['g'] = ref
	params['o'] = outfile

	cmd = '{{args.snvsniffer}} somatic %s "%s" "%s" "{{i.tumor}}" "{{i.normal}}"' % (cmdargs(params), theader, nheader)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## strelka
	{% when 'strelka' %}
	# config
	configParams = {{args.configParams}}

	configParams['normalBam']      = {{i.normal | quote}}
	configParams['tumorBam']       = {{i.tumor | quote}}
	configParams['referenceFasta'] = ref
	configParams['runDir']         = {{job.outdir | quote}}

	cmd = '{{args.strelka}} %s' % cmdargs(configParams)
	runcmd (cmd)

	# run
	params['m'] = 'local'
	params['g'] = mem2({{args.mem | quote}}, 'G')[:-1]
	params['j'] = {{args.nthread}}
	cmd = '{{job.outdir}}/runWorkflow.py %s' % cmdargs(params)
	runcmd (cmd)
	# merge
	mem = mem2 ({{args.mem | quote}}, 'java')
	mergeparams = {
		'R'                   : ref,
		'V:SNV'               : "{{job.outdir}}/results/variants/somatic.snvs.vcf.gz",
		'V:INDEL'             : "{{job.outdir}}/results/variants/somatic.indels.vcf.gz",
		'o'                   : outfile,
		'genotypeMergeOptions': 'UNIQUIFY'
	}

	cmd = '{{args.gatk}} -T CombineVariants %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(mergeparams, dash = '-', equal = ' '))
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## virmid
	{% when 'virmid' %}
	mem = mem2 ({{args.mem | quote}}, 'java')
	params['R'] = ref
	params['D'] = {{i.tumor | quote}}
	params['N'] = {{i.normal | quote}}
	params['w'] = {{job.outdir | quote}}
	cmd = '{{args.virmid}} %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params))
	runcmd (cmd)
	runcmd ('mv "{{job.outdir}}/*.virmid.som.passed.vcf" "%s"' % outfile)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## vardict
	{% when 'vardict' %}
	params['v'] = True
	params['G'] = ref
	params['b'] = "{{i.tumor}}|{{i.normal}}"

	cmd = '{{args.vardict}} %s > "%s"' % (cmdargs(params), outfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

{% endcase %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
