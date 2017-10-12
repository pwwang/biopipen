from os import path
from shutil import move
from sys import stderr
from collections import OrderedDict

ref  = {{args.ref | quote}}
	
{{runcmd}}
{{params2CmdArgs}}

params = OrderedDict({{args.params}})
fq1 = {{out.fq1 | quote}}
fq2 = {{out.fq2 | quote}}
try:
	{% if args.tool | lambda x: x == 'wgsim' %}
	{% if args.gz %}
	fq1 = "{{out.fq1 | [:-3]}}"
	fq2 = "{{out.fq2 | [:-3]}}"
	{% endif %}

	params['N'] = {{args.num}}
	params['1'] = {{args.len1}}
	params['2'] = {{args.len2}}
	params['S'] = {{args.seed | lambda x: -1 if x is None else x}}
	cmd = '{{args.wgsim}} %s "%s" "%s" "%s"' % (params2CmdArgs(params, noq = ['N', '1', '2', 'S']), ref, fq1, fq2)
	runcmd (cmd)

	{% if args.gz %}
	runcmd ('gzip "%s"' % fq1)
	runcmd ('gzip "%s"' % fq2)
	{% endif %}
	
	{% elif args.tool | lambda x: x == 'dwgsim' %}
	prefix = {{out.fq1 | [:-8] | quote}}
	{% if args.gz %}
	fq1 = "{{out.fq1 | [:-3]}}"
	fq2 = "{{out.fq2 | [:-3]}}"
	prefix = "{{out.fq1 | [:-11]}}"
	{% endif %}

	params['N'] = {{args.num}}
	params['1'] = {{args.len1}}
	params['2'] = {{args.len2}}
	{% if args.seed | lambda x: x is not None %}
	params['z'] = {{args.seed}}
	{% endif %}
	params['S'] = {{args.seed | lambda x: 0 if x is None else x}}
	cmd = '{{args.dwgsim}} %s "%s" "%s"' % (params2CmdArgs(params, noq = ['N', '1', '2', 'S']), ref, prefix)
	runcmd (cmd)

	if path.exists (prefix + '.bwa.read1.fastq.gz'):
		move (prefix + '.bwa.read1.fastq.gz', prefix + '_1.fastq.gz')
		move (prefix + '.bwa.read2.fastq.gz', prefix + '_2.fastq.gz')
		{% if args.gz | lambda x: not x %}
		runcmd ('gunzip "%s"' % (prefix + '_1.fastq.gz'))
		runcmd ('gunzip "%s"' % (prefix + '_2.fastq.gz'))
		{% endif %}
	else:
		move (prefix + '.bwa.read1.fastq', prefix + '_1.fastq')
		move (prefix + '.bwa.read2.fastq', prefix + '_2.fastq')
		{% if args.gz %}
		runcmd ('gzip "%s"' % (prefix + '_1.fastq'))
		runcmd ('gzip "%s"' % (prefix + '_2.fastq'))
		{% endif %}
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise