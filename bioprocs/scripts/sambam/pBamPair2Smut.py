from os import path, makedirs
from shutil import rmtree
from sys import stderr

if not path.exists ("{{bring.tumor[0]}}"):
	stderr.write ("Input file '{{in._tumor}}' is not indexed.")
	exit (1)
if not path.exists ("{{bring.normal[0]}}"):
	stderr.write ("Input file '{{in._normal}}' is not indexed.")
	exit (1)

if not path.exists ({{args.ref | quote}}):
	stderr.write ("Reference file not specified")
	exit (1)

{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{in.tumor | fn | fn}}.{{in.normal | fn | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

ref     = {{args.ref | quote}}
tool    = {{args.tool | quote}}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{out.outfile | quote}}
params  = {{args.params}}
if gz:	outfile = outfile[:-3]
try:
	############## gatk
	{% if args.tool | lambda x: x == 'gatk' %}
	
	intvfile = "{{job.outdir}}/interval.list"
	cmd      = '{{args.samtools}} idxstats "{{in.tumor}}" | head -n -1 | cut -f1 > "%s"' % (intvfile)
	runcmd (cmd)

	mem = mem2 ({{args.mem | quote}}, 'java')

	params['R']        = ref
	params['I:tumor']  = {{in.tumor | quote}}
	params['I:normal'] = {{in.normal | quote}}
	params['o']        = outfile
	params['nct']      = {{args.nthread}}
	params['L']        = intvfile

	cmd = '{{args.gatk}} -T MuTect2 %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, dash = '-', equal = ' '))
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## somaticsniper
	{% elif args.tool | lambda x: x == 'somaticsniper' %}
	params['f'] = ref
	params['F'] = 'vcf'
	
	cmd = '{{args.somaticsniper}} %s "{{in.tumor}}" "{{in.normal}}" "%s"' % (params2CmdArgs(params), outfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))
	
	############## snvsniffer
	{% elif args.tool | lambda x: x == 'snvsniffer' %}
	theader = "{{job.outdir}}/{{in.tumor | bn}}.header"
	cmd = '{{args.samtools}} view -H "{{in.tumor}}" > "%s"' % theader
	runcmd (cmd)
	nheader = "{{job.outdir}}/{{in.normal | bn}}.header"
	cmd = '{{args.samtools}} view -H "{{in.normal}}" > "%s"' % nheader
	runcmd (cmd)

	params['g'] = ref
	params['o'] = outfile

	cmd = '{{args.snvsniffer}} somatic %s "%s" "%s" "{{in.tumor}}" "{{in.normal}}"' % (params2CmdArgs(params), theader, nheader)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))
	
	############## strelka
	{% elif args.tool | lambda x: x == 'strelka' %}
	# config
	configParams = {{args.configParams}}

	configParams['normalBam']      = {{in.normal | quote}}
	configParams['tumorBam']       = {{in.tumor | quote}}
	configParams['referenceFasta'] = ref
	configParams['runDir']         = {{job.outdir | quote}}

	cmd = '{{args.strelka}} %s' % params2CmdArgs(configParams)
	runcmd (cmd)
	
	# run
	params['m'] = 'local'
	params['g'] = mem2({{args.mem | quote}}, 'G')[:-1]
	params['j'] = {{args.nthread}}
	cmd = '{{job.outdir}}/runWorkflow.py %s' % params2CmdArgs(params)
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

	cmd = '{{args.gatk}} -T CombineVariants %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(mergeparams, dash = '-', equal = ' '))
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## virmid
	{% elif args.tool | lambda x: x == 'virmid' %}
	mem = mem2 ({{args.mem | quote}}, 'java')
	params['R'] = ref
	params['D'] = {{in.tumor | quote}}
	params['N'] = {{in.normal | quote}}
	params['w'] = {{job.outdir | quote}}
	cmd = '{{args.virmid}} %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params))
	runcmd (cmd)
	runcmd ('mv "{{job.outdir}}/*.virmid.som.passed.vcf" "%s"' % outfile)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############## vardict
	{% elif args.tool | lambda x: x == 'vardict' %}
	params['v'] = True
	params['G'] = ref
	params['b'] = "{{in.tumor}}|{{in.normal}}"
	
	cmd = '{{args.vardict}} %s > "%s"' % (params2CmdArgs(params), outfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)