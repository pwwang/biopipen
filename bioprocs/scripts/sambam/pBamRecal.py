
from os import path, makedirs, remove, symlink
from sys import stderr
from shutil import rmtree, move
from time import sleep

if not path.exists ({{bring.infile[0] | quote}}):
	stderr.write ("Input file '{{in._infile}}' is not indexed.")
	exit (1)

if not path.exists ({{args.knownSites | quote}}) and {{args.tool | quote}} == 'gatk':
	stderr.write ("knownSites file is required by GATK but is not specified (args.knownSites) or not exists.")
	exit (1)
	
{{ runcmd }}
{{ mem2 }}
{{ buildrefIndex }}
{{ params2CmdArgs }}

ref = {{args.ref | quote}}

tmpdir    = {{ args.tmpdir | quote }}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

params = {{args.params}}
try:
	########## gatk
	{% if args.tool | lambda x: x == 'gatk' %}
	mem                          = mem2({{ args.mem | quote }}, 'java')
	intfile                      = "{{out.outfile | prefix}}.intervals"
	paramsRealignerTargetCreator = {{args.paramsRealignerTargetCreator}}

	paramsRealignerTargetCreator['R'] = ref
	paramsRealignerTargetCreator['I'] = {{in.infile | quote}}
	paramsRealignerTargetCreator['o'] = intfile

	runcmd ('{{args.gatk}} -T RealignerTargetCreator %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(paramsRealignerTargetCreator, dash = '-', equal = ' ')))

	bamfileIr            = "{{out.outfile | prefix}}.ir.bam"
	paramsIndelRealigner = {{args.paramsIndelRealigner}}

	paramsIndelRealigner['R']               = ref
	paramsIndelRealigner['I']               = {{in.infile | quote}}
	paramsIndelRealigner['o']               = bamfileIr
	paramsIndelRealigner['targetIntervals'] = intfile

	runcmd ('{{args.gatk}} -T IndelRealigner %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(paramsIndelRealigner, dash = '-', equal = ' ')))

	recaltable             = "{{out.outfile | prefix}}.recaltable"
	paramsBaseRecalibrator = {{args.paramsBaseRecalibrator}}

	paramsBaseRecalibrator['R']          = ref
	paramsBaseRecalibrator['I']          = bamfileIr
	paramsBaseRecalibrator['o']          = recaltable
	paramsBaseRecalibrator['knownSites'] = {{args.knownSites}}

	runcmd ('{{args.gatk}} -T BaseRecalibrator %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(paramsBaseRecalibrator, dash = '-', equal = ' ')))

	paramsPrintReads = {{args.paramsPrintReads}}

	paramsPrintReads['R']    = ref
	paramsPrintReads['I']    = bamfileIr
	paramsPrintReads['o']    = {{out.outfile | quote}}
	paramsPrintReads['BQSR'] = recaltable

	runcmd ('{{args.gatk}} -T PrintReads %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, params2CmdArgs(paramsPrintReads, dash = '-', equal = ' ')))
	remove (bamfileIr)
	move ("{{out.outfile | prefix}}.bai", "{{out.idxfile}}")
	
	########## bamutil
	{% elif args.tool | lambda x: x == 'bamutil' %}
	{% if args.knownSites %}
	params['dbsnp'] = {{args.knownSites | quote}}
	{% endif %}
	params['in']      = {{in.infile | quote}}
	params['out']     = {{out.outfile | quote}}
	params['refFile'] = ref

	cmd     = '{{args.bamutil}} recab %s'
	ref2    = "{{ proc.workdir }}/0/{{ args.ref | bn }}"
	params2 = {k:v for k,v in params.items()}

	params2['refFile'] = ref2

	cmd1 = cmd % params2CmdArgs(params, equal = ' ')
	cmd2 = cmd % params2CmdArgs(params2, equal = ' ')

	refcache = "{{ args.ref | prefix }}-bs.umfa"
	if path.exists (refcache):
		runcmd (cmd1)
	else:
		r = buildrefIndex (
			{{job.index}},
			ref,
			[refcache],
			cmd1,
			{{proc.workdir | quote}},
			cmd2
		)
		params2['refFile'] = r
		{% if job.index | lambda x: x > 0 %}
		runcmd(cmd % params2CmdArgs(params2, equal = ' '))
		{% endif %}
	
	cmd = '{{args.samtools}} index "{{out.outfile}}" "{{out.idxfile}}"'
	runcmd (cmd)
	
	{% endif %}

except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)