from sys import stderr

	
{{runcmd}}
{{params2CmdArgs}}
params = {{args.params}}

qcdir = {{in.qcdir | quote}}
try:
	{% if args.tool | lambda x: x == 'multiqc' %}
	params['o'] = {{out.outdir | quote}}
	params['p'] = True
	cmd = '{{args.multiqc}} "{{in.qcdir}}" %s' % params2CmdArgs(params)
	runcmd(cmd)

	{% else %}
	raise Exception('Tool {{args.tool}} not supported.')

	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise