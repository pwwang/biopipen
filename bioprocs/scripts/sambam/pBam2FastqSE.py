from os import makedirs, path

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

infile  = {{in.infile | quote}}
fqfile  = {{out.fqfile | quote}}
{% if args.gz %}
fqfile = fqfile[:-3]
{% endif %}

{{runcmd}}
{{mem2}}
{{params2CmdArgs}}

params  = {{args.params}}
try:
	############# biobambam
	{% if args.tool | lambda x: x == 'biobambam' %}
	params['gz']       = 0
	#bug
	#params['S']        = fqfile
	params['filename'] = infile
	params['T']        = path.join(tmpdir, infile + '.tmp')
	if infile.endswith('.sam'):
		params['inputformat'] = 'sam'
	
	cmd = '{{args.biobambam}} %s > "%s"' % (params2CmdArgs(params, dash = '', equal = '='), fqfile)
	runcmd (cmd)

	############# bedtools
	{% elif args.tool | lambda x: x == 'bedtools' %}
	params['i']  = infile
	params['fq'] = fqfile

	cmd = '{{args.bedtools}} bamtofastq %s' % params2CmdArgs(params, dash = '-', equal = ' ')
	runcmd (cmd)

	############# samtools
	{% elif args.tool | lambda x: x == 'samtools' %}
	params['t'] = True
	params['s'] = fqfile

	cmd = '{{args.samtools}} fastq %s "%s"' % (params2CmdArgs(params), infile)
	runcmd (cmd)

	############# picard
	{% elif args.tool | lambda x: x == 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'java')
	params['TMP_DIR'] = tmpdir
	params['I'] = infile
	params['F'] = fqfile
	cmd = '{{args.picard}} SamToFastq %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(params, dash = '', equal = '='))
	runcmd (cmd)
	{% endif %}

	{% if args.gz %}
	runcmd ('gzip "%s"' % (fqfile))
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)