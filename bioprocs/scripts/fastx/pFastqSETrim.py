from sys import stderr
from shutil import move

from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs

params = {{args.params}}
try:
{% case args.tool %}
	{% when 'trimmomatic' %}
	mem    = mem2 ({{args.mem | quote}}, "java")
	minlen = str({{args.minlen}} * 2)
	adfile = "{{job.outdir}}/adapters.fa"
	with open (adfile, "w") as ad:
		ad.write (">TruSeq3_IndexedAdapter\n")
		ad.write ({{args.adapter | quote}} + "\n")

	params['threads'] = {{args.nthread}}
	cmd    = '{{args.trimmomatic}} %s SE %s "{{in.fq}}" "{{out.outfq}}" ILLUMINACLIP:%s:2:30:10 LEADING:{{args.cut5}} TRAILING:{{args.cut3}} SLIDINGWINDOW:4:{{args.minq}} MINLEN:%s' % (mem, cmdargs(params, dash = '-', equal = ' '), adfile, minlen)
	runcmd (cmd)

	{% when 'cutadapt' %}
	params['a'] = {{args.adapter | quote}}
	params['u'] = "{{args.cut5}}"
	params['u'] = "-{{args.cut3}}"
	params['m'] = {{args.minlen}}
	params['q'] = "{{args.minq}},{{args.minq}}"
	params['o'] = {{out.outfq | quote}}
	cmd = '{{args.cutadapt}} %s "{{in.fq}}"' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)

	{% when 'skewer' %}
	params['m'] = 'any'
	params['t'] = {{args.nthread}}
	params['x'] = {{args.adapter | quote}}
	params['Q'] = {{args.minq}}
	params['l'] = {{args.minlen}}
	params['z'] = {{args.gz}}
	params['o'] = "{{job.outdir}}/tmp"
	cmd = '{{args.skewer}} %s "{{in.fq}}"' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)
	outfq = "{{job.outdir}}/tmp-trimmed.fastq"
	{% if args.gz %}
	outfq += ".gz"
	{% endif %}
	move (outfq, "{{out.outfq}}")

	{% else %}
	raise Exception ('Tool {{args.tool}} not supported')
{% endcase %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
