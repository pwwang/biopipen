from bioprocs.utils import logger, runcmd, cmdargs
from pyppl import Box

params = {{args.params}}

qcdir = {{in.qcdir | quote}}
try:
	{% if args.tool | lambda x: x == 'multiqc' %}
	params['o'] = {{out.outdir | quote}}
	params['p'] = True
	cmd = '{{args.multiqc}} "{{in.qcdir}}" %s' % cmdargs(params)
	runcmd(cmd)

	{% else %}
	raise Exception('Tool {{args.tool}} not supported.')

	{% endif %}
except Exception as ex:
	logger.error("Job failed: %s" % str(ex))
	raise
