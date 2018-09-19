import shlex
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.parallel import Parallel

# vcf2maf doesn't support multiple samples
# so we have to split multi-sample vcf to single-sample vcfs
# convert each single-sample vcf to maf and then combine them
{% if args.tool == 'vcf2maf' %}
from os import path
from subprocess import check_output
vep = check_output(['which', {{args.vep | quote}}]).strip()

params = {{args.params}}

# get samples
{% if args.samfunc %}
samples = ({{args.samfunc}})({{i.infile | quote}})
{% else %}
samples = check_output([{{args.bcftools | quote}}, 'query', '-l', {{i.infile | quote}}]).splitlines()
samples = list(filter(None, [x.strip() for x in samples]))
{% endif %}

{% 	if args.somatic %}
params['input-vcf']  = {{i.infile | quote}}
params['output-maf'] = {{o.outfile | quote}}
params['vep-data']   = {{args.vepDb | quote}}
params['vep-forks']  = {{args.nthread}}
params['filter-vcf'] = {{args.filtervcf | quote}}
params['ref-fasta']  = {{args.ref | quote}}
params['vep-path']   = path.dirname(vep)
{% if args.tumor1st %}
params['tumor-id']   = samples.pop(0)
params['normal-id']  = samples[0] if samples else 'NORMAL'
{% else %}
params['normal-id']  = samples.pop(0)
params['tumor-id']   = samples[0] if samples else 'NORMAL'
{% endif %}

cmd = '{{args.vcf2maf}} %s' % (cmdargs(params, equal=' '))
runcmd(cmd)

{% 	else %}
cmds = []
for sample in samples:
	vtparams = {}
	vtparams['a'] = True
	vtparams['c'] = sample
	vtparams['e'] = True
	samplevcf     = "{{job.outdir}}/{{i.infile | fn}}-%s.vcf" % sample
	cmd = '{{args.vcftools}} %s {{i.infile | quote}} > "%s"' % (cmdargs(vtparams), samplevcf)

	# vcf2maf.pl --input-vcf ZYYP-ZYYB.vcf  --output-maf ZYYP-ZYYB.snpEff.maf --tumor-id ZXLT-ZXLB_TUMOR --normal-id ZXLT-ZXLB_NORMAL --vep-data /path/to/vep/cache/ --filter-vcf /path/to/vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --ref-fasta /path/to/hs37d5/phase2_reference_assembly_sequence/hs37d5.fa --vep-path /path/to/miniconda2/bin
	params['input-vcf']  = samplevcf
	params['output-maf'] = "{{job.outdir}}/{{i.infile | fn}}-%s.maf" % sample
	params['vep-data']   = {{args.vepDb | quote}}
	params['vep-forks']  = {{args.nthread}}
	params['filter-vcf'] = {{args.filtervcf | quote}}
	params['ref-fasta']  = {{args.ref | quote}}
	params['vep-path']   = path.dirname(vep)

	cmd = cmd + '; {{args.vcf2maf}} --tumor-id %s %s' % (sample, cmdargs(params, equal=' '))
	cmds.append(cmd)

{% 		if args.nthread == 1 %}
for cmd in cmds: runcmd(cmd)
{% 		else %}
# Note the threads may be hanging on here.
p = Parallel({{args.nthread}})
p.run('{}', [(cmd,) for cmd in cmds])
{% 		endif %}

for i, sample in enumerate(samples):
	singlemaf = "{{job.outdir}}/{{i.infile | fn}}-%s.maf" % sample
	if i == 0:
		runcmd('cat "%s" > {{o.outfile | quote}}' % singlemaf)
	else:
		runcmd('egrep -v "^#|^Hugo_Symbol" "%s" >> {{o.outfile | quote}}' % singlemaf)
{% 	endif %}
{% endif %}