from os import path
from glob import glob

{{runcmd}}
{{params2CmdArgs}}

invcfs = glob(path.join({{in.indir | quote}}, {{args.pattern | quote}}))
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
cmds.append(cmd)

{% endif %}