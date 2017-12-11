import shlex
{{ runcmd }}
{{ params2CmdArgs }}
{{ parallel }}

# vcf2maf doesn't support multiple samples
# so we have to split multi-sample vcf to single-sample vcfs
# convert each single-sample vcf to maf and then combine them
{% if args.tool | lambda x: x == 'vcf2maf' %}
from os import path
from subprocess import check_output
vep = check_output(['which', {{args.vep | quote}}]).strip()

params = {{args.params}}

# get samples
samples = check_output([{{args.bcftools | quote}}, 'query', '-l', {{in.infile | quote}}]).splitlines()
samples = list(filter(None, [x.strip() for x in samples]))

{% 	if args.somatic %}
params['input-vcf']  = {{in.infile | quote}}
params['output-maf'] = {{out.outfile | quote}}
params['vep-data']   = {{args.vepDb | quote}}
params['filter-vcf'] = {{args.filtervcf | quote}}
params['ref-fasta']  = {{args.ref | quote}}
params['vep-path']   = path.dirname(vep)
params['tumor-id']   = samples.pop(0)
params['normal-id']  = samples[0] if samples else 'NORMAL'

cmd = '{{args.vcf2maf}} %s' % (params2CmdArgs(params, equal=' '))
runcmd(cmd)

{% 	else %}
cmds = []
for sample in samples:
	vtparams = {}
	vtparams['a'] = True
	vtparams['c'] = sample
	vtparams['e'] = True
	samplevcf     = "{{job.outdir}}/{{in.infile | fn}}-%s.vcf" % sample
	cmd = '{{args.vcftools}} %s {{in.infile | quote}} > "%s"' % (params2CmdArgs(vtparams), samplevcf)

	# vcf2maf.pl --input-vcf ZYYP-ZYYB.vcf  --output-maf ZYYP-ZYYB.snpEff.maf --tumor-id ZXLT-ZXLB_TUMOR --normal-id ZXLT-ZXLB_NORMAL --vep-data /data2/junwenwang/shared/reference/hg19/vep/cache/ --filter-vcf /data2/junwenwang/shared/reference/hg19/vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --ref-fasta /data2/junwenwang/shared/reference/hs37d5/phase2_reference_assembly_sequence/hs37d5.fa --vep-path /data2/junwenwang/shared/tools/miniconda2/bin
	params['input-vcf']  = samplevcf
	params['output-maf'] = "{{job.outdir}}/{{in.infile | fn}}-%s.maf" % sample
	params['vep-data']   = {{args.vepDb | quote}}
	params['filter-vcf'] = {{args.filtervcf | quote}}
	params['ref-fasta']  = {{args.ref | quote}}
	params['vep-path']   = path.dirname(vep)

	cmd = cmd + '; {{args.vcf2maf}} --tumor-id %s %s' % (sample, params2CmdArgs(params, equal=' '))
	cmds.append(cmd)

{% 		if args.nthread | lambda x: x == 1 %}
for cmd in cmds: runcmd(cmd)
{% 		else %}
parallel(cmds, {{args.nthread}})
{% 		endif %}

for i, sample in enumerate(samples):
	singlemaf = "{{job.outdir}}/{{in.infile | fn}}-%s.maf" % sample
	if i == 0:
		runcmd('cat "%s" > {{out.outfile | quote}}' % singlemaf)
	else:
		runcmd('egrep -v "^#|^Hugo_Symbol" "%s" >> {{out.outfile | quote}}' % singlemaf)
{% 	endif %}
{% endif %}