from os import path
from sys import stderr
from glob import glob

{{runcmd}}
{{params2CmdArgs}}
{{parallel}}

invcfs = {{in.infiles}}

cmds     = []
invcfgzs = []
for invcf in invcfs:
	if invcf.endswith('.gz'): 
		invcfgzs.append(invcf)
		continue
	invcfgz = path.join('{{job.outdir}}', path.basename(invcf) + '.gz')
	invcfgzs.append(invcfgz)
	cmd = 'bgzip "%s" -c > "%s"; tabix -p vcf "%s"' % (invcf, invcfgz, invcfgz)
	cmds.append(cmd)
invcfs = invcfgzs

{% if args.nthread | lambda x: x > 1 %}
parallel(cmds, {{args.nthread}})
{% else %}
for cmd in cmds: runcmd(cmd)
{% endif %}

########### vcftools
{% if args.tool | lambda x: x == 'vcftools' %}
params = {}
params['d'] = True
params['t'] = True
cmd = '{{args.vcftools}} %s %s > "{{out.outfile}}"' % (params2CmdArgs(params, equal=' '), ' '.join(invcfs))
runcmd(cmd)

########### gatk
{% elif args.tool | lambda x: x == 'gatk' %}
params = {}
params['o'] = {{out.outfile | quote}}
params['R'] = {{args.ref | quote}}
for i, invcf in enumerate(invcfs):
	key = 'variant' + ' '*i
	params[key] = invcf
params['genotypemergeoption'] = 'UNIQUIFY'

cmd = '{{args.gatk}} -T CombineVariants %s' % (params2CmdArgs(params, equal=' '))
runcmd(cmd)

{% endif %}