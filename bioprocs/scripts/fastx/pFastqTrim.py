from sys import stderr
from shutil import move

from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs

params = {{args.params}}
try:
	{% if args.tool | lambda x: x == 'trimmomatic' %}
	def seqrev (seq):
		d = {
			'A':'T', 'T':'A', 'G':'C', 'C':'G',
			'a':'t', 't':'a', 'g':'c', 'c':'g'
		}
		return ''.join([d[s] for s in seq])

	mem    = mem2 ({{args.mem | quote}}, 'java')
	minlen = str({{args.minlen}} * 2)
	adfile = "{{job.outdir}}/adapters.fa"
	with open (adfile, "w") as ad:
		ad.write (">PE1\n")
		ad.write (seqrev({{args.adapter1 | quote}}) + "\n")
		ad.write (">PE1_rc\n")
		ad.write ({{args.adapter1 | quote}} + "\n")
		ad.write (">PE2\n")
		ad.write (seqrev({{args.adapter2 | quote}}) + "\n")
		ad.write (">PE2_rc\n")
		ad.write ({{args.adapter2 | quote}} + "\n")

	params['threads'] = {{args.nthread}}
	cmd = '{{args.trimmomatic}} %s PE %s "{{in.fq1}}" "{{in.fq2}}" "{{out.outfq1}}" /dev/null "{{out.outfq2}}" /dev/null ILLUMINACLIP:%s:2:30:10 LEADING:{{args.cut5}} TRAILING:{{args.cut3}} SLIDINGWINDOW:4:{{args.minq}} MINLEN:%s' % (mem, cmdargs(params, dash = '-', equal = ' '), adfile, minlen)
	runcmd (cmd)

	{% elif args.tool | lambda x: x == 'cutadapt' %}
	params['a'] = {{args.adapter1 | quote}}
	params['A'] = {{args.adapter2 | quote}}
	params['u'] = "{{args.cut5}}"
	params['u'] = "-{{args.cut3}}"
	params['U'] = "{{args.cut5}}"
	params['U'] = "-{{args.cut3}}"
	params['m'] = {{args.minlen}}
	params['q'] = "{{args.minq}},{{args.minq}}"
	params['o'] = {{out.outfq1 | quote}}
	params['p'] = {{out.outfq2 | quote}}
	cmd = '{{args.cutadapt}} %s {{ in.fq1 | quote }} {{ in.fq2 | quote }}' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)

	{% elif args.tool | lambda x: x == 'skewer' %}
	params['m'] = 'pe'
	params['t'] = {{args.nthread}}
	params['x'] = {{args.adapter1 | quote}}
	params['y'] = {{args.adapter1 | quote}}
	params['Q'] = {{args.minq}}
	params['l'] = {{args.minlen}}
	params['z'] = {{args.gz}}
	params['o'] = "{{job.outdir}}/tmp"
	cmd = '{{args.skewer}} %s "{{in.fq1}}" "{{in.fq2}}"' % cmdargs(params, dash = '-', equal = ' ')
	runcmd (cmd)
	outfq1 = "{{job.outdir}}/tmp-trimmed-pair1.fastq"
	outfq2 = "{{job.outdir}}/tmp-trimmed-pair2.fastq"
	{% if args.gz %}
	outfq1 += ".gz"
	outfq2 += ".gz"
	{% endif %}
	move (outfq1, "{{out.outfq1}}")
	move (outfq2, "{{out.outfq2}}")

	{% else %}
	raise Exception ('{{args.tool}} not supported')
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
