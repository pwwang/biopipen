import re
from os import path, symlink
from sys import stdout, stderr
from time import sleep
from shutil import move
from collections import OrderedDict

{{ runcmd }}
{{ buildrefIndex }}
{{ params2CmdArgs }}
		
# determine whether the reference file is specified
ref = "{{args.ref}}"
if not ref or not path.exists (ref):
	stderr.write ('No reference specified or found.')
	exit (1)

# detemine default read group
rg = {{ args.rg }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{out.outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{out.outfile | fn}}.L{{job.index}}"
if not rg['SM']:
	rg['SM'] = "{{out.outfile | fn}}"

params = OrderedDict({{args.params}})

def sam2bam(samfile, bamfile):
	cmd = '{{args.samtools}} view -Sb "%s" > "%s"; rm -f "%s"' % (samfile, bamfile, samfile)
	runcmd(cmd)

############ bowtie2
{% if args.tool | lambda x: x == 'bowtie2' %}
# check reference index files
ref = buildrefIndex (
	{{job.index}}, 
	ref, 
	["{{ args.ref }}.1.bt2", "{{ args.ref }}.2.bt2", "{{ args.ref }}.3.bt2", "{{ args.ref }}.4.bt2", "{{ args.ref }}.rev.1.bt2", "{{ args.ref }}.rev.2.bt2"], 
	'{{args.bowtie2_build}} --thread {{args.nthread}} "{{args.ref}}" "{{args.ref}}"',
	"{{ proc.workdir }}",
	'{{args.bowtie2_build}} --thread {{args.nthread}} "{{ job.outdir }}/{{ args.ref | bn }}" "{{ job.outdir }}/{{ args.ref | bn }}"')

# do mapping
params['threads'] = {{args.nthread}}
params['rg-id'] = rg['ID']
i = 0
for k,v in rg.items():
	if k == 'ID': continue
	params['rg' + ' '*i] = k + ':' + v
	i += 1
params['x'] = ref
params['1'] = {{in.fq1 | quote}}
params['2'] = {{in.fq2 | quote}}
params['S'] = {% if args.outfmt | lambda x: x == 'bam' %}{{out.outfile | prefix | lambda x: x + '.sam' | quote}}{% else %}{{out.outfile | quote}}{% endif %} 

cmd = '{{args.bowtie2}} %s' % params2CmdArgs(params, equal = ' ', noq = ['threads'])
runcmd (cmd)

{% if args.outfmt | lambda x: x == 'bam' %}
sam2bam(params['S'], {{out.outfile | quote}})
{% endif %}

############ bwa
{% elif args.tool | lambda x: x == 'bwa' %}
# check reference index
ref = buildrefIndex (
	{{job.index}}, 
	ref, 
	["{{ args.ref }}.bwt", "{{ args.ref }}.amb", "{{ args.ref }}.ann", "{{ args.ref }}.pac", "{{ args.ref }}.sa"], 
	'{{args.bwa}} index "{{args.ref}}"',
	"{{ proc.workdir }}",
	'{{args.bwa}} index "{{ job.outdir }}/{{ args.ref | bn }}"')
		
# do mapping
params['t'] = {{args.nthread}}
params['R'] = "@RG\\tID:%s\\t%s" % (rg['ID'], "\\t".join(k + ':' +v for k,v in rg.items() if k!='ID'))
outfile     = {% if args.outfmt | lambda x: x == 'bam' %}{{out.outfile | prefix | lambda x: x + '.sam' | quote}}{% else %}{{out.outfile | quote}}{% endif %} 
cmd = '{{args.bwa}} mem %s "%s" {{in.fq1 | quote}} {{in.fq2 | quote}} > "%s"' % (params2CmdArgs(params, noq=['t']), ref, outfile)
runcmd (cmd)

{% if args.outfmt | lambda x: x == 'bam' %}
sam2bam(outfile, {{out.outfile | quote}})
{% endif %}

{% elif args.tool | lambda x: x == 'ngm' %}
# check reference index
ref = buildrefIndex (
	{{job.index}}, 
	ref, 
	["{{ args.ref }}-enc.2.ngm", "{{ args.ref }}-ht-13-2.3.ngm"], 
	'{{args.ngm}} -r "{{args.ref}}" -t {{args.nthread}}',
	"{{ proc.workdir }}",
	'{{args.ngm}} -r "{{ job.outdir }}/{{ args.ref | bn }}" -t {{args.nthread}}')
	
# do mapping
params['1'] = {{in.fq1 | quote}}
params['2'] = {{in.fq2 | quote}}
params['r'] = ref
params['o'] = {{out.outfile | quote}}
params['b'] = {{ args.outfmt | lambda x: x=='bam' }}
params.update({('rg-' + k.lower()):v for k,v in rg.items()})
params['t'] = {{args.nthread}}
cmd = '{{args.ngm}} %s' % (params2CmdArgs(params, equal = ' ', noq = ['t']))
runcmd (cmd)

############ star
{% elif args.tool | lambda x: x=='star' %}
# check reference index
ref = buildrefIndex (
	{{job.index}}, 
	ref, 
	map(lambda x: path.join("{{args.ref | prefix}}.star", x), ["chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt", "exonGeTrInfo.tab", "exonInfo.tab", "geneInfo.tab", "Genome", "genomeParameters.txt", "SA", "SAindex", "sjdbInfo.txt", "sjdbList.fromGTF.out.tab", "sjdbList.out.tab", "transcriptInfo.tab"]), 
	'mkdir -p "{{args.ref | prefix}}.star"; {{args.star}} --runMode genomeGenerate --genomeDir "{{args.ref | prefix}}.star" --genomeFastaFiles "{{args.ref}}" --sjdbOverhang 100 --sjdbGTFfile "{{args.refgene}}" --runThreadN {{args.nthread}}',
	"{{ proc.workdir }}",
	'mkdir -p "{{job.outdir}}/{{args.ref | fn}}.star"; {{args.star}} --runMode genomeGenerate --genomeDir "{{job.outdir}}/{{args.ref | fn}}.star" --genomeFastaFiles "{{job.outdir}}/{{args.ref | bn}}" --sjdbOverhang 100 --sjdbGTFfile "{{args.refgene}}" --runThreadN {{args.nthread}}')
refDir = "{{args.ref | prefix}}.star" if ref == "{{args.ref}}" else "{{job.outdir}}/{{args.ref | fn}}.star"
rfcmd  = "zcat" if "{{in.fq1}}".endswith('.gz') else 'bzcat' if "{{in.fq1}}".endswith('.bz2') else 'cat'

params['genomeDir'] = refDir
params['readFilesIn'] = '{{in.fq1 | quote}} {{in.fq2 | quote}}'
params['readFilesCommand'] = rfcmd
params['readNameSeparator'] = '.'
params['outFileNamePrefix'] = "{{job.outdir}}/"
params['outSAMtype'] = {{args.outfmt.upper() | quote}}
cmd = '{{args.star}} %s' % params2CmdArgs(params, equal = ' ', noq = ['readFilesIn', 'readNameSeparator'])
runcmd (cmd)

outfile = "{{job.outdir}}/Aligned.out.{{args.outfmt}}"
if path.exists(outfile):
	move (outfile, "{{out.outfile}}")
{% endif %}