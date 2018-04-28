
from os import path, makedirs, remove, symlink
from sys import stderr
from shutil import rmtree, move
from time import sleep

from pyppl import Box
from bioprocs.utils import runcmd, mem2, cmdargs, logger, reference
from bioprocs.utils.poll import Poll

# check index
'''
if not path.exists ({ {bring.infile[0] | quote}}):
	stderr.write ("Input file '{ {in._infile}}' is not indexed.")
	exit (1)
'''

if not path.exists ({{args.knownSites | quote}}) and {{args.tool | quote}} == 'gatk':
	stderr.write ("knownSites file is required by GATK but is not specified (args.knownSites) or not exists.")
	exit (1)

workdir  = {{proc.workdir | quote}}
outdir   = {{job.outdir | quote}}
ref      = {{args.ref | quote}}
ref2     = path.join(outdir, path.basename(ref))
joblen   = {{proc.size}}
jobindex = {{job.index}}
tmpdir   = {{ args.tmpdir | quote }}
tmpdir   = path.join (tmpdir, "{{proc.id}}.{{in.infile | fn}}.{{job.index}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

poll = Poll(workdir, joblen, jobindex)

def checkAndBuildIndex(indexes, buildcmd, buildcmd2):
	def todo(indexes, ref, buildcmd, ref2, buildcmd2):
		logger.info('Checking if reference file is indexed.')
		if not reference.checkIndex(indexes):
			logger.info('No, trying to build index at the first job.')
			reference.buildIndex(ref, buildcmd, ref2, buildcmd2)
	poll.first(todo, indexes, ref, buildcmd, ref2, buildcmd2)
	return ref2 if path.exists(ref2) else ref

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

	runcmd ('{{args.gatk}} -T RealignerTargetCreator %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(paramsRealignerTargetCreator, dash = '-', equal = ' ')))

	bamfileIr            = "{{out.outfile | prefix}}.ir.bam"
	paramsIndelRealigner = {{args.paramsIndelRealigner}}

	paramsIndelRealigner['R']               = ref
	paramsIndelRealigner['I']               = {{in.infile | quote}}
	paramsIndelRealigner['o']               = bamfileIr
	paramsIndelRealigner['targetIntervals'] = intfile

	runcmd ('{{args.gatk}} -T IndelRealigner %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(paramsIndelRealigner, dash = '-', equal = ' ')))

	recaltable             = "{{out.outfile | prefix}}.recaltable"
	paramsBaseRecalibrator = {{args.paramsBaseRecalibrator}}

	paramsBaseRecalibrator['R']          = ref
	paramsBaseRecalibrator['I']          = bamfileIr
	paramsBaseRecalibrator['o']          = recaltable
	paramsBaseRecalibrator['knownSites'] = {{args.knownSites}}

	runcmd ('{{args.gatk}} -T BaseRecalibrator %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(paramsBaseRecalibrator, dash = '-', equal = ' ')))

	paramsPrintReads = {{args.paramsPrintReads}}

	paramsPrintReads['R']    = ref
	paramsPrintReads['I']    = bamfileIr
	paramsPrintReads['o']    = {{out.outfile | quote}}
	paramsPrintReads['BQSR'] = recaltable

	runcmd ('{{args.gatk}} -T PrintReads %s -Djava.io.tmpdir="%s" %s' % (mem, tmpdir, cmdargs(paramsPrintReads, dash = '-', equal = ' ')))
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

	cmd1 = cmd % cmdargs(params, equal = ' ')
	cmd2 = cmd % cmdargs(params2, equal = ' ')
	refcache = "{{ args.ref | prefix }}-bs.umfa"

	if path.exists(refcache):
		runcmd(cmd1)
	else:
		r = checkAndBuildIndex([refcache], cmd1, cmd2)

		params2['refFile'] = r
		{% if job.index | lambda x: x > 0 %}
		runcmd(cmd % cmdargs(params2, equal = ' '))
		{% endif %}

	cmd = '{{args.samtools}} index "{{out.outfile}}" "{{out.idxfile}}"'
	runcmd (cmd)

	{% endif %}

except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
