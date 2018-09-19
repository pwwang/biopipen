from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stdout, stderr
from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs

infile    = {{ i.infile | quote }}
outfile   = {{ o.outfile | quote }}
tmpdir    = {{ args.tmpdir | quote }}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{i.infile | fn}}.{{job.index}}")
doRmdup   = {{ args.rmdup }}

if not path.exists (tmpdir):
	makedirs (tmpdir)

params = {{args.params}}
try:
{% case args.tool %}
	########### biobambam
	{% when 'biobambam'%}
	if doRmdup:
		params['rmdup'] = 1
		params['D']     = "/dev/null"
	params['I'] = infile
	params['O'] = outfile
	params['tmpfile'] = path.join(tmpdir, 'tmp.')
	params['markthreads'] = {{args.nthread}}

	cmd = '{{args.biobambam_bamsort}} %s' % cmdargs(params, dash='', equal = '=')
	runcmd (cmd)

	########### sambamba
	{% when 'sambamba' %}
	if doRmdup:
		params['r'] = True
	params['t'] = {{args.nthread}}
	params['tmpdir'] = tmpdir

	cmd = '{{args.sambamba}} markdup %s "%s" "%s"' % (cmdargs(params), infile, outfile)
	runcmd (cmd)

	########### samtools
	{% when 'samtools' %}
	cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, outfile)
	runcmd (cmd)

	########### picard
	{% when 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')

	params['-Djava.io.tmpdir'] = tmpdir
	if doRmdup:
		params['REMOVE_DUPLICATES'] = "true"
	params['TMP_DIR'] = tmpdir
	params['I'] = infile
	params['O'] = outfile
	params['M'] = "/dev/null"

	mfile = "/dev/null"
	cmd = '{{args.picard}} MarkDuplicates %s %s %s ' % (mem, cmdargs(params, dash='', equal='='), rmdup)
	runcmd (cmd)

	########### bamutil
	{% when 'bamutil' %}
	if doRmdup:
		params['rmDups'] = True
	params['in']  = infile
	params['out'] = outfile
	cmd = '{{args.bamutil}} dedup %s' % cmdargs(params, equal = ' ')
	runcmd (cmd)

{% endcase %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
