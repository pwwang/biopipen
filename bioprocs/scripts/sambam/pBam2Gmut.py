from os import path, makedirs
from shutil import rmtree
from sys import stderr, exit
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, mem2

tmpdir    = path.join ("{{ args.tmpdir}}", "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

ref     = {{args.ref | quote}}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{o.outfile | quote}}
if gz:	outfile = outfile[:-3]

params = {{args.params}}
try:
{% case args.tool %}
	############ gatk
	{% when 'gatk' %}
	mem = mem2({{args.mem | quote}}, 'java')

	params['R']   = ref
	params['I']   = {{i.infile | quote}}
	params['o']   = outfile
	params['nct'] = {{args.nthread}}

	cmd = '{{args.gatk}} -T HaplotypeCaller %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '-', equal = ' '))
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############ vardict
	{% when 'vardict' %}
	params['v'] = True
	params['G'] = ref
	params['b'] = {{i.infile | quote}}

	cmd = '{{args.vardict}} %s > "%s"' % (cmdargs(params), outfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############ snvsniffer
	{% when 'snvsniffer' %}
	hfile = "{{job.outdir}}/{{i.infile | bn}}.header"
	cmd   = '{{args.samtools}} view -H "{{i.infile}}" > "%s"' % hfile
	runcmd (cmd)

	params['g'] = ref
	params['o'] = outfile

	cmd = '{{args.snvsniffer}} snp %s "%s" "{{i.infile}}"' % (cmdargs(params), hfile)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############ platypus
	{% when 'platypus' %}
	params['refFile']     = ref
	params['bamFiles']    = {{i.infile | quote}}
	params['nCPU']        = {{args.nthread}}
	params['output']      = outfile
	params['logFileName'] = outfile + '.log'

	cmd = '{{args.platypus}} callVariants %s' % cmdargs(params)
	runcmd (cmd)
	if gz:	runcmd ('gzip "%s"' % (outfile))

	############ strelka
	{% when 'strelka' %}
	# config
	configParams = {{args.configParams}}
	configParams['bam']            = {{i.infile | quote}}
	configParams['referenceFasta'] = ref
	configParams['runDir']         = {{job.outdir | quote}}

	cmd = '{{args.strelka}} %s' % cmdargs(configParams)
	runcmd (cmd)

	# run
	params['m'] = 'local'
	params['j'] = {{args.nthread}}
	params['g'] = mem2({{args.mem | quote}}, 'G')[:-1]

	cmd = '{{job.outdir}}/runWorkflow.py %s' % cmdargs(params)
	runcmd (cmd)
	ofile = "{{job.outdir}}/results/variants/genome.S1.vcf.gz"
	runcmd ('mv "%s" "%s.gz"' % (ofile, outfile))
	if not gz: runcmd ('gunzip "%s.gz"' % outfile)

{% endcase %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
