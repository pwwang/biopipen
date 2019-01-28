from os import makedirs, path
from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs, logger

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

# bam2fastq will create {i.infile}.tmp, use file in indir in case of permission issue
infile  = {{i.infile | quote}}
fqfile1 = {{o.fqfile1 | quote}}
fqfile2 = {{o.fqfile2 | quote}}
{% if args.gz %}
fqfile1 = fqfile1[:-3]
fqfile2 = fqfile2[:-3]
{% endif %}

params  = {{args.params}}
try:
{% case args.tool %}
{% when 'biobambam' %}
	params['gz'] = 0
	params['F']  = fqfile1
	params['F2'] = fqfile2
	params['T']  = path.join(tmpdir, infile + '.tmp')
	params['filename'] = infile
	if infile.endswith('.sam'):
		params['inputformat'] = 'sam'
	cmd = '{{args.biobambam}} %s' % cmdargs(params, dash = '', equal = '=')
	runcmd (cmd)
{% when 'bedtools' %}
	params['i']   = infile
	params['fq']  = fqfile1
	params['fq2'] = fqfile2
	cmd = '{{args.bedtools}} bamtofastq %s' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)
{% when 'samtools' %}
	params['t'] = True
	params['1'] = fqfile1
	params['2'] = fqfile2
	cmd = '{{args.samtools}} fastq %s "%s"' % (cmdargs(params, dash = '-', equal = ' '), infile)
	runcmd (cmd)
{% when 'picard' %}
	mem = mem2({{ args.mem | quote }}, 'Java')
	params[mem]                = True
	params['-Djava.io.tmpdir'] = tmpdir
	params['TMP_DIR']          = tmpdir
	params['I']                = infile
	params['F']                = fqfile1
	params['F2']               = fqfile2
	cmd = '{{args.picard}} SamToFastq %s' % cmdargs(params, dash='', equal='=')
	runcmd (cmd)
{% endcase %}

{% if args.gz %}
	runcmd ('gzip "%s"' % (fqfile1))
	runcmd ('gzip "%s"' % (fqfile2))
{% endif %}
except Exception as ex:
	logger.error ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)
