from sys import stderr
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

fq   = {{i.fq | quote}}
params = {{args.params}}
try:
	{% if args.tool == 'fastqc' %}
	params['o'] = {{o.outdir | quote}}
	cmd = '{{args.fastqc}} %s "{{i.fq}}"' % cmdargs(params)
	runcmd(cmd)

	{% else %}
	raise Exception('Tool {{args.tool}} %s not supported.')

	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise