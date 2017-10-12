from os import makedirs, path

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

infile  = {{in.infile | quote}}
fqfile1 = {{out.fqfile1 | quote}}
fqfile2 = {{out.fqfile2 | quote}}
{% if args.gz %}
fqfile1 = fqfile1[:-3]
fqfile2 = fqfile2[:-3]
{% endif %}

{{runcmd}}
{{mem2}}
{{params2CmdArgs}}

params  = {{args.params}}
try:
	{% if args.tool | lambda x: x == 'biobambam' %}
		params['gz'] = 0
		params['F']  = fqfile1
		params['F2'] = fqfile2
		params['T']  = path.join(tmpdir, infile + '.tmp')
		params['filename'] = infile
		if infile.endswith('.sam'):
			params['inputformat'] = 'sam'
		cmd = '{{args.biobambam_bamtofastq}} %s' % params2CmdArgs(params, dash = '', equal = '=', noq = 'gz')
		runcmd (cmd)
	{% elif args.tool | lambda x: x == 'bedtools' %}
		params['i']   = infile
		params['fq']  = fqfile1
		params['fq2'] = fqfile2
		cmd = '{{args.bedtools}} bamtofastq %s' % params2CmdArgs(params, dash = '-', equal = ' ')
		runcmd (cmd)
	{% elif args.tool | lambda x: x == 'samtools' %}
		params['t'] = True
		params['1'] = fqfile1
		params['2'] = fqfile2
		cmd = '{{args.samtools}} fastq %s "%s"' % (params2CmdArgs(params, dash = '-', equal = ' ', noq = ['t']), infile)
		runcmd (cmd)
	{% elif args.tool | lambda x: x == 'picard' %}
		mem = mem2({{ args.mem | quote }}, 'Java')
		params[mem]                = True
		params['-Djava.io.tmpdir'] = tmpdir
		params['TMP_DIR']          = tmpdir
		params['I']                = infile
		params['F']                = fqfile1
		params['F2']               = fqfile2
		cmd = '{{args.picard}} SamToFastq %s' % params2CmdArgs(params, dash='', equal='=')
		runcmd (cmd)
	{% endif %}
		
	{% if args.gz %}
	runcmd ('gzip "%s"' % (fqfile1))
	runcmd ('gzip "%s"' % (fqfile2))
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)