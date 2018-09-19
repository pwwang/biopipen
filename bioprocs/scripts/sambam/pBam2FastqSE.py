from os import makedirs, path
from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

infile  = {{i.infile | quote}}
fqfile  = {{o.fqfile | quote}}
{% if args.gz %}
fqfile = fqfile[:-3]
{% endif %}

params  = {{args.params}}
try:
{% case args.tool %}
	############# biobambam
	{% when 'biobambam' %}
	params['gz']       = 0
	#bug
	#params['S']        = fqfile
	params['filename'] = infile
	params['T']        = path.join(tmpdir, infile + '.tmp')
	if infile.endswith('.sam'):
		params['inputformat'] = 'sam'

	cmd = '{{args.biobambam}} %s > "%s"' % (cmdargs(params, dash = '', equal = '='), fqfile)
	runcmd (cmd)

	############# bedtools
	{% when 'bedtools' %}
	params['i']  = infile
	params['fq'] = fqfile

	cmd = '{{args.bedtools}} bamtofastq %s' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)

	############# samtools
	{% when 'samtools' %}
	params['t'] = True
	params['s'] = fqfile

	cmd = '{{args.samtools}} fastq %s "%s"' % (cmdargs(params), infile)
	runcmd (cmd)

	############# picard
	{% when 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	params['TMP_DIR'] = tmpdir
	params['I'] = infile
	params['F'] = fqfile
	cmd = '{{args.picard}} SamToFastq %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '', equal = '='))
	runcmd (cmd)
{% endcase %}

	{% if args.gz %}
	runcmd ('gzip "%s"' % (fqfile))
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)
