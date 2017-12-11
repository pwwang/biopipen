from sys import stderr
from subprocess import check_output
{{runcmd}}
{{params2CmdArgs}}

allsamples = check_output([{{args.bcftools | quote}}, 'query', '-l', {{in.infile | quote}}]).splitlines()
allsamples = list(filter(None, [x.strip() for x in allsamples]))

samples = {{in.samples | lambda x: list(filter(None, [x.strip() for x in x.split(',')]))}}
if not samples: samples = allsamples

nonexists = [sample for sample in samples if sample not in allsamples]
if nonexists:
	stderr.write('pyppl.log.warning: Samples not exist: %s\n' % nonexists)
	for ne in nonexists: del samples[samples.index(ne)]

cmds = []
########### vcftools
{% if args.tool | lambda x: x == 'vcftools' %}
for sample in samples:
	params = {}
	params['a'] = True
	params['c'] = sample
	params['e'] = True
	cmd = '{{args.vcftools}} %s "{{in.infile}}" > "{{out.outdir}}/%s.vcf"' % (params2CmdArgs(params), sample)
	cmds.append(cmd)

########### gatk
{% elif args.tool | lambda x: x == 'gatk' %}
for sample in samples:
	params                       = {}
	params['R']                  = {{args.ref | quote}}
	params['V']                  = {{in.infile | quote}}
	params['o']                  = "{{out.outdir}}/%s.vcf" % sample
	params['sample_name']        = sample
	params['excludeFiltered']    = True
	params['excludeNonVariants'] = True

	cmd = '{{args.gatk}} -T SelectVariants %s' % (params2CmdArgs(params, equal=' '))
	cmds.append(cmd)

{% endif %}

{% if args.nthread | lambda x: x == 1 %}
for cmd in cmds: runcmd(cmd)
{% else %}
{{parallel}}
parallel(cmds, {{args.nthread}})
{% endif %}
