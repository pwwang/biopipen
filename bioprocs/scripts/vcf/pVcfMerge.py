from os import path
from sys import stderr
from glob import glob

from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.parallel import Parallel

invcfs  = {{i.infiles | repr}}
nthread = int({{args.nthread | repr}})

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

if nthread > 1:
	#parallel(cmds, {{args.nthread}})
	p = Parallel(nthread, raiseExc = True)
	p.run('{}', [(cmd, ) for cmd in cmds])
else:
	for cmd in cmds: runcmd(cmd)

########### vcftools
{% if args.tool == 'vcftools' %}
params = {}
params['d'] = True
params['t'] = True
cmd = '{{args.vcftools}} %s %s > "{{o.outfile}}"' % (cmdargs(params, equal=' '), ' '.join(invcfs))
runcmd(cmd)

########### gatk
{% elif args.tool == 'gatk' %}
params = {}
params['o'] = {{o.outfile | quote}}
params['R'] = {{args.ref | quote}}
for i, invcf in enumerate(invcfs):
	key = 'variant' + ' '*i
	params[key] = invcf
params['genotypemergeoption'] = 'UNIQUIFY'

cmd = '{{args.gatk}} -T CombineVariants %s' % (cmdargs(params, equal=' '))
runcmd(cmd)

{% endif %}