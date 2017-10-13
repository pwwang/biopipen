import re
from os import path, symlink, remove
from sys import stdout, stderr, exit
from time import sleep
from collections import OrderedDict

{{ runcmd }}
{{ buildrefIndex }}
{{ params2CmdArgs }}
		
# determine whether the reference file is specified
ref = "{{args.ref}}"

# detemine default read group
rg = {{ args.rg }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{out.outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{out.outfile | fn}}.L{{job.index}}"
if not rg['SM']:
	rg['SM'] = "{{out.outfile | fn}}"

params = OrderedDict({{args.params}})

############# bowtie2
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
	params['rg' + ' '*i] = k.lower() + ':' + v
	i += 1
params['x'] = ref
params['U'] = {{in.fq | quote}}
params['S'] = {{out.outfile | quote}}

cmd = '{{args.bowtie2}} %s' % params2CmdArgs(params, equal = ' ', noq = ['threads'])
runcmd (cmd)

############# bwa
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

cmd = '{{args.bwa}} mem %s "%s" "{{in.fq}}" > "{{out.outfile}}"' % (params2CmdArgs(params, noq=['t']), ref)
runcmd (cmd)

############# ngm
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
params['q'] = {{in.fq | quote}}
params['r'] = ref
params['o'] = {{out.outfile | quote}}
params['b'] = {{ args.outformat | lambda x: x=='bam' }}
params.update({('rg-' + k.lower()):v for k,v in rg.items()})
params['t'] = {{args.nthread}}

cmd = '{{args.ngm}} %s' % (params2CmdArgs(params, equal = ' ', noq = ['t']))
runcmd (cmd)
{% endif %}
		