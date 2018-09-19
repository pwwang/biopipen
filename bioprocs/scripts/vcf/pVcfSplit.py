from sys import stderr
from pyppl import Box
from subprocess import check_output
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.parallel import Parallel

allsamples = check_output([{{args.bcftools | quote}}, 'query', '-l', {{i.infile | quote}}]).splitlines()
allsamples = list(filter(None, [x.strip() for x in allsamples]))

samples = {{i.samples | lambda x: list(filter(None, [x.strip() for x in x.split(',')]))}}
if not samples: samples = allsamples

nonexists = [sample for sample in samples if sample not in allsamples]
if nonexists:
	stderr.write('pyppl.log.warning: Samples not exist: %s\n' % nonexists)
	for ne in nonexists: del samples[samples.index(ne)]

cmds = []
########### vcftools
{% if args.tool == 'vcftools' %}
for sample in samples:
	params = {}
	params['a'] = True
	params['c'] = sample
	params['e'] = True
	params.update({{args.params}})
	cmd = '{{args.vcftools}} %s "{{i.infile}}" > "{{o.outdir}}/{{i.infile | fn}}-%s.vcf"' % (cmdargs(params), sample)
	cmds.append(cmd)

########### awk
{% elif args.tool == 'awk' %}
# write the awk script
awkfile = '{{job.outdir}}/vcfsample.awk'
awkfh   = open(awkfile, 'w')
awkfh.write("""
BEGIN {
	OFS="\\t"
} 
$0 ~ "^##" {
	print
} 
$0 ~ "^#CHROM" {
	print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t"sample
}	
$0 !~ "^#" {
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,$index
}
""")
awkfh.close()
for i, sample in enumerate(samples):
	cmd = '{{args.awk}} -v sample="{sample}" index={index} -f {awkfile} {infile}'.format(
		sample = sample,
		index  = 10 + i,
		awk    = str(repr(awkfile)),
		infile = str(repr(infile))
	)
	cmds.append(cmd)

########### gatk
{% elif args.tool == 'gatk' %}
for sample in samples:
	params                       = {}
	params['R']                  = {{args.ref | quote}}
	params['V']                  = {{i.infile | quote}}
	params['o']                  = "{{o.outdir}}/{{i.infile | fn}}-%s.vcf" % sample
	params['sample_name']        = sample
	params['excludeFiltered']    = True
	params['excludeNonVariants'] = True
	params.update({{args.params}})
	cmd = '{{args.gatk}} -T SelectVariants %s' % (cmdargs(params, equal=' '))
	cmds.append(cmd)

{% endif %}

{% if args.nthread == 1 %}
for cmd in cmds: runcmd(cmd)
{% else %}
p = Parallel({{args.nthread}})
p.run('{}', [(cmd,) for cmd in cmds])
{% endif %}
