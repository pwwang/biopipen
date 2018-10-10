
from os import makedirs, path
from shutil import rmtree
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2

tmpdir    = path.join ({{ args.tmpdir | quote}}, "{{proc.id}}.{{i.infiles[0] | fn2}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

params = {{args.params}}
try:
{% case args.tool %}
	############# picard
	{% when 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	infiles = {{ i.infiles }}
	for i, infile in enumerate(infiles):
		params['I' + ' ' * i] = infile
	{% if args.nthread > 1 %}
	params['USE_THREADING'] = 'true'
	{% else %}
	params['USE_THREADING'] = 'false'
	{% endif %}
	params['TMP_DIR'] = tmpdir
	params['O']       = {{o.outfile | quote}}
	params['AS']      = 'true'

	cmd = '{{args.picard}} MergeSamFiles %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '', equal = '='))
	runcmd (cmd)

	############# bamutil
	{% when 'bamutil' %}
	infiles = {{ i.infiles }}
	for i, infile in enumerate(infiles):
		params['i' + ' ' * i] = infile
	params['o'] = {{o.outfile | quote}}

	cmd = '{{args.bamutil}} mergeBam %s' % cmdargs(params)
	runcmd (cmd)

	############# samtools
	{% when 'samtools' %}
	inlist = path.join({{job.outdir | quote}}, 'bamlist.txt')
	with open(inlist, 'w') as f:
		f.write('\n'.join({{i.infiles}}) + '\n')
	params['@'] = {{args.nthread}}
	params['O'] = 'bam'
	params['b'] = inlist

	cmd = '{{args.samtools}} merge %s {{o.outfile | quote}}' % cmdargs(params)
	runcmd (cmd)

	############# sambamba
	{% when 'sambamba' %}
	params['t'] = {{args.nthread}}
	cmd = '{{args.sambamba}} merge %s {{o.outfile | quote}} {{ i.infiles | asquote }}' % cmdargs(params)
	runcmd (cmd)

{% endcase %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
