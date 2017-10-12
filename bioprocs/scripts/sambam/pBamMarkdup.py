from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stdout, stderr

{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}

infile    = {{ in.infile | quote }}
outfile   = {{ out.outfile | quote }}
tmpdir    = {{ args.tmpdir | quote }}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
doRmdup   = {{ args.rmdup }}

if not path.exists (tmpdir):
	makedirs (tmpdir)

params = {{args.params}}
try:
	########### biobambam
	{% if args.tool | lambda x: x == 'biobambam'%}
	if doRmdup:
		params['rmdup'] = 1
		params['D']     = "/dev/null"
	params['I'] = infile
	params['O'] = outfile
	params['tmpfile'] = path.join(tmpdir, 'tmp.')
	params['markthreads'] = {{args.nthread}}

	cmd = '{{args.biobambam_bamsort}} %s' % params2CmdArgs(params, dash='', equal = '=', noq = ['rmdup', 'D', 'markthreads'])
	runcmd (cmd)

	########### sambamba
	{% elif args.tool | lambda x: x == 'sambamba' %}
	if doRmdup:
		params['r'] = True
	params['t'] = {{args.nthread}}
	params['tmpdir'] = tmpdir

	cmd = '{{args.sambamba}} markdup %s "%s" "%s"' % (params2CmdArgs(params, noq=['t']), infile, outfile)
	runcmd (cmd)
	
	########### samtools
	{% elif args.tool | lambda x: x == 'samtools' %}
	cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, outfile)
	runcmd (cmd)

	########### picard
	{% elif args.tool | lambda x: x == 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	
	params['-Djava.io.tmpdir'] = tmpdir
	if doRmdup:
		params['REMOVE_DUPLICATES'] = "true"
	params['TMP_DIR'] = tmpdir
	params['I'] = infile
	params['O'] = outfile
	params['M'] = "/dev/null"

	mfile = "/dev/null"
	cmd = '{{args.picard}} MarkDuplicates %s %s %s ' % (mem, params2CmdArgs(params, dash='', equal='=', noq=['M', 'REMOVE_DUPLICATES']), rmdup)
	runcmd (cmd)

	########### bamutil
	{% elif args.tool | lambda x: x == 'bamutil' %}
		if doRmdup:
			params['rmDups'] = True
		params['in']  = infile
		params['out'] = outfile
		cmd = '{{args.bamutil}} dedup %s' % params2CmdArgs(params, equal = ' ')
		runcmd (cmd)
	
	{% endif %}
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)