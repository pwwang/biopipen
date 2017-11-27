
from os import makedirs, path
from shutil import rmtree

tmpdir    = path.join ({{ args.tmpdir | quote}}, "{{proc.id}}.{{in.infiles | fsDirname}}.{{job.index}}")
if not path.exists (tmpdir): makedirs (tmpdir)
	
{{ runcmd }}
{{ mem2 }}
{{ params2CmdArgs }}

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

	cmd = '{{args.picard}} MergeSamFiles %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, dash = '', equal = '='))
	runcmd (cmd)

	############# bamutil
	{% elif args.tool | lambda x: x == 'bamutil' %}
	infiles = {{ in.infiles }}
	for i, infile in enumerate(infiles):
		params['i' + ' ' * i] = infile
	params['o'] = {{out.outfile | quote}}

	cmd = '{{args.bamutil}} mergeBam %s' % params2CmdArgs(params)
	runcmd (cmd)
	
	############# samtools
	{% elif args.tool | lambda x: x == 'samtools' %}
	inlist = path.join({{job.outdir | quote}}, 'bamlist.txt')
	with open(inlist, 'w') as f:
		f.write('\n'.join({{in.infiles}}) + '\n')
	params['@'] = {{args.nthread}}
	params['O'] = 'bam'
	params['b'] = inlist

	cmd = '{{args.samtools}} merge %s {{out.outfile | quote}}' % params2CmdArgs(params)
	runcmd (cmd)

	############# sambamba
	{% elif args.tool | lambda x: x == 'sambamba' %}
	params['t'] = {{args.nthread}}
	cmd = '{{args.sambamba}} merge %s {{out.outfile | quote}} {{ in.infiles | asquote }}' % params2CmdArgs(params)
	runcmd (cmd)

	{% endif %}
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)