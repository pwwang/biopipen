
from os import makedirs, path
from shutil import rmtree
from pyppl import Box
from bioprocs.utils import cmdargs, runcmd, mem2

tmpdir    = path.join ({{ args.tmpdir | quote}}, "{{proc.id}}.{{in.infiles.0 | fn}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)

params = {{args.params}}
try:
	############# picard
	{% if args.tool | lambda x: x == 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	infiles = {{ in.infiles }}
	for i, infile in enumerate(infiles):
		params['I' + ' ' * i] = infile
	{% if args.nthread | lambda x: x>1 %}
	params['USE_THREADING'] = 'true'
	{% else %}
	params['USE_THREADING'] = 'false'
	{% endif %}
	params['TMP_DIR'] = tmpdir
	params['O']       = {{out.outfile | quote}}
	params['AS']      = 'true'

	cmd = '{{args.picard}} MergeSamFiles %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(params, dash = '', equal = '='))
	runcmd (cmd)

	############# bamutil
	{% elif args.tool | lambda x: x == 'bamutil' %}
	infiles = {{ in.infiles }}
	for i, infile in enumerate(infiles):
		params['i' + ' ' * i] = infile
	params['o'] = {{out.outfile | quote}}

	cmd = '{{args.bamutil}} mergeBam %s' % cmdargs(params)
	runcmd (cmd)

	############# samtools
	{% elif args.tool | lambda x: x == 'samtools' %}
	inlist = path.join({{job.outdir | quote}}, 'bamlist.txt')
	with open(inlist, 'w') as f:
		f.write('\n'.join({{in.infiles}}) + '\n')
	params['@'] = {{args.nthread}}
	params['O'] = 'bam'
	params['b'] = inlist

	cmd = '{{args.samtools}} merge %s {{out.outfile | quote}}' % cmdargs(params)
	runcmd (cmd)

	############# sambamba
	{% elif args.tool | lambda x: x == 'sambamba' %}
	params['t'] = {{args.nthread}}
	cmd = '{{args.sambamba}} merge %s {{out.outfile | quote}} {{ in.infiles | asquote }}' % cmdargs(params)
	runcmd (cmd)

	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
