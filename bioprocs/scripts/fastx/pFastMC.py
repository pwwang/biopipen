from bioprocs.utils import logger, runcmd, cmdargs
from pyppl import Box

params = {{args.params}}

qcdir = {{i.qcdir | quote}}
try:
	{% if args.tool == 'multiqc' %}
	params['o'] = {{o.outdir | quote}}
	params['p'] = True
	cmd = '{{args.multiqc}} "{{i.qcdir}}" %s' % cmdargs(params)
	runcmd(cmd)

	{% else %}
	raise Exception('Tool {{args.tool}} not supported.')

	{% endif %}
except Exception as ex:
	logger.error("Job failed: %s" % str(ex))
	raise
