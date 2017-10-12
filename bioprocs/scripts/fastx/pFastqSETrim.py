from sys import stderr
from shutil import move

{{runcmd}}
{{mem2}}
{{params2CmdArgs}}

params = {{args.params}}
try:
	{% if args.tool | lambda x: x=='trimmomatic' %}
	mem    = mem2 ({{args.mem | quote}}, "java")
	minlen = str({{args.minlen}} * 2)
	adfile = "{{job.outdir}}/adapters.fa"
	with open (adfile, "w") as ad:
		ad.write (">TruSeq3_IndexedAdapter\n")
		ad.write ({{args.adapter | quote}} + "\n")

	params['threads'] = {{args.nthread}}
	cmd    = '{{args.trimmomatic}} %s SE %s "{{in.fq}}" "{{out.outfq}}" ILLUMINACLIP:%s:2:30:10 LEADING:{{args.cut5}} TRAILING:{{args.cut3}} SLIDINGWINDOW:4:{{args.minq}} MINLEN:%s' % (mem, params2CmdArgs(params, dash = '-', equal = ' ', noq = ['threads']), adfile, minlen)
	runcmd (cmd)

	{% elif args.tool | lambda x: x=='cutadapt' %}
	params['a'] = {{args.adapter | quote}}
	params['u'] = "{{args.cut5}}"
	params['u'] = "-{{args.cut3}}"
	params['m'] = {{args.minlen}}
	params['q'] = "{{args.minq}},{{args.minq}}"
	params['o'] = {{out.outfq | quote}}
	cmd = '{{args.cutadapt}} %s "{{in.fq}}"' % params2CmdArgs(params, dash = '-', equal = ' ', noq = ['-u', '-m'])
	runcmd (cmd)

	{% elif args.tool | lambda x: x=='skewer' %}
	params['m'] = 'any'
	params['t'] = {{args.nthread}}
	params['x'] = {{args.adapter | quote}}
	params['Q'] = {{args.minq}}
	params['l'] = {{args.minlen}}
	params['z'] = {{args.gz}}
	params['o'] = "{{job.outdir}}/tmp"
	cmd = '{{args.skewer}} %s "{{in.fq}}"' % params2CmdArgs(params, dash = '-', equal = ' ', noq = ['m', 't', 'x', 'Q', 'l'])
	runcmd (cmd)
	outfq = "{{job.outdir}}/tmp-trimmed.fastq"
	{% if args.gz %}
	outfq += ".gz"
	{% endif %}
	move (outfq, "{{out.outfq}}")
			
	{% else %}
	raise Exception ('Tool {{args.tool}} not supported')
	{% endif %}
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise