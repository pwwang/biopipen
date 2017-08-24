"""
A set of processes to generate/process sam/bam files
"""

from pyppl import proc
from .utils import mem, runcmd, buildArgIndex, checkArgsRef, checkArgsRefgene, buildArgsFastaFai, buildArgsFastaDict, polling0, pollingAll

"""
@name:
	pSam2Bam
@description:
	Deal with mapped sam/bam files, including sort, markdup, and/or index
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output bam file
	`idxfile:file`: The index of the output bam file
	- If args.index == False, it'll a link to outfile and should be never used
@args:
	`tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools)
	`sambamba`         : The path of the sambamba. Default: sambamba 
	`picard`           : The path of the picard. Default: picard 
	`biobambam_bamsort`: The path of the biobambam's bamsort. Default: bamsort 
	`samtools`         : The path of the samtools. Default: samtools 
	`sort`             : Do sorting? Default: True 
	- If input is sam, tool is biobambam, this should be True
	`index`            : Do indexing? Default: True
	`markdup`          : Do duplicates marking? Default: False
	- `rmdup` for samtools will be called
	`rmdup`            : Do duplicates removing? Default: False
	`tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
	`sortby`           : Sort by coordinate or queryname. Default: coordinate
	`nthread`          : Default: 1
	`informat`         : The format of input file. Default: <detect from extension> (sam|bam)
	`params`           : Other parameters for `tool`. Defaut: ""
	`mem`              : The max memory to use. Default: "16G"
	- Unit could be G/g/M/m
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
@requires:
	[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
"""
pSam2Bam                        = proc (desc = 'Deal with mapped sam/bam files, including sort, markdup, rmdup, and/or index.')
pSam2Bam.input                  = "infile:file"
pSam2Bam.output                 = "outfile:file:{{infile | fn}}.bam, idxfile:file:{{infile | fn}}.bam.bai"
pSam2Bam.args.tool              = "biobambam" 
pSam2Bam.args.sambamba          = "sambamba"
pSam2Bam.args.picard            = "picard"
pSam2Bam.args.biobambam_bamsort = "bamsort"
pSam2Bam.args.samtools          = "samtools"
pSam2Bam.args.sort              = True
pSam2Bam.args.index             = True
pSam2Bam.args.markdup           = False
pSam2Bam.args.rmdup             = False
pSam2Bam.args.tmpdir            = __import__('tempfile').gettempdir() 
pSam2Bam.args.sortby            = "coordinate" 
pSam2Bam.args.nthread           = 1
pSam2Bam.args.informat          = ""
pSam2Bam.args.params            = ""
pSam2Bam.args.mem               = "16G"	
pSam2Bam.args._mem2M            = mem.toM.python
pSam2Bam.args._memtoJava        = mem.toJava.python
pSam2Bam.args._runcmd           = runcmd.python
pSam2Bam.lang                   = "python"
pSam2Bam.script                 = """
from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stdout, stderr, exit

{{ args._runcmd }}
{{ args._mem2M }}
{{ args._memtoJava }}

infile    = {{ infile | quote }}
outfile   = {{ outfile | quote }}
tool      = {{ args.tool | quote }}
informat  = {{ args.informat | quote }}
informat  = informat if informat else {{ infile | ext | [1:] | quote}} 
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{infile | fn}}.{{#}}")
doSort    = {{ args.sort }}
doIndex   = {{ args.index }}
doMarkdup = {{ args.markdup }}
doRmdup   = {{ args.rmdup }}
if doRmdup:
	doMarkdup = True
if not path.exists (tmpdir):
	makedirs (tmpdir)
try:
	if tool == 'biobambam':
		mem = memtoM({{ args.mem | quote }})
		indexArg = ""
		if doIndex:
			indexArg = 'index=1 indexfilename="{{idxfile}}"'
		cmd = '{{args.biobambam_bamsort}} I="%s" O="%s" SO={{args.sortby}} blockme="%s" tmpfile="%s/tmp." inputformat="%s" outformat=bam inputthreads={{args.nthread}} outputthreads={{args.nthread}} %s markduplicates=%s rmdup=%s {{args.params}}' % (infile, outfile, mem, tmpdir, informat, indexArg, str(int(doMarkdup)), str(int(doRmdup)))
		runcmd (cmd)
	elif tool == 'sambamba':
		if not (doSort or doIndex or doMarkdup or doRmdup):
			cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (outfile, infile)
			runcmd (cmd)
		else:
			bamfile = outfile
			if informat == 'sam':
				bamfile = "{{job.outdir}}/{{infile | fn}}.s2b.bam"
				cmd = '{{args.sambamba}} view -S -f bam -o %s -t {{args.nthread}} %s' % (bamfile, infile)
				runcmd (cmd)
				infile = bamfile
			if doSort:
				sortby = ''
				if {{args.sortby | quote}} == 'queryname':
					sortby = '-n -N'
				bamfile = "{{job.outdir}}/{{infile | fn}}.sorted.bam"
				cmd = '{{args.sambamba}} sort -m {{args.mem}} --tmpdir="%s" -o "%s" %s -t {{args.nthread}} {{args.params}} "%s"' % (tmpdir, bamfile, sortby, infile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doMarkdup:
				rmdup = ""
				if doRmdup:
					rmdup = "-r"
				bamfile = "{{job.outdir}}/{{infile | fn}}.dedup.bam"
				cmd = '{{args.sambamba}} markdup %s -t {{args.nthread}} --tmpdir="%s" "%s" "%s"' % (rmdup, tmpdir, infile, bamfile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doIndex:
				if path.exists (infile + '.bai'):
					move (infile + '.bai', {{idxfile | quote}})
				else:
					cmd = '{{args.sambamba}} index -t {{args.nthread}} "%s" "%s"' % (infile, {{idxfile | quote}})
					runcmd (cmd)
			if infile != outfile:
				if path.exists(infile + '.bai'):
					move (infile + '.bai', outfile + '.bai')
				move (infile, outfile)
	elif tool == 'samtools':
		if not (doSort or doIndex or doMarkdup or doRmdup):
			cmd = '{{args.samtools}} view -b -o "%s" -O bam "%s"' % (outfile, infile)
			runcmd (cmd)
		else:
			bamfile = outfile
			if doSort:
				mem = memtoM({{ args.mem | quote }})
				sortby = ''
				if {{args.sortby | quote}} == 'queryname':
					sortby = '-n'
				bamfile = "{{job.outdir}}/{{infile | fn}}.sorted.bam" 
				cmd = '{{args.samtools}} sort -m %sM %s -o "%s" -T "%s" -@ {{args.nthread}} -O bam "%s"' % (mem, sortby, bamfile, tmpdir, infile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doMarkdup or doRmdup:
				bamfile = "{{job.outdir}}/{{infile | fn}}.dedup.bam"
				cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, bamfile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doIndex:
				cmd = '{{args.samtools}} index "%s" "%s"' % (bamfile, {{idxfile | quote}})
				runcmd (cmd)
			if infile != outfile:
				if path.exists(infile + '.bai'):
					move (infile + '.bai', outfile + '.bai')
				move (infile, outfile)
	elif tool == 'picard':
		mem = memtoJava({{ args.mem | quote }})
		if not (doSort or doIndex or doMarkdup or doRmdup):
			cmd = '{{args.picard}} SamFormatConverter %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, outfile)
			runcmd (cmd)
		else:
			bamfile = outfile
			if doSort:
				bamfile = "{{job.outdir}}/{{infile | fn}}.sorted.bam" 
				cmd = '{{args.picard}} SortSam %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" SO={{args.sortby}}' % (mem, tmpdir, tmpdir, infile, bamfile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doMarkdup:
				rmdup = ""
				if doRmdup:
					rmdup = "REMOVE_DUPLICATES=true"
				mfile = "/dev/null"
				bamfile = "{{job.outdir}}/{{infile | fn}}.dedup.bam"
				cmd = '{{args.picard}} MarkDuplicates %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" M="%s" ' % (mem, tmpdir, tmpdir, infile, bamfile, mfile)
				runcmd (cmd)
				if infile != {{infile | quote}}:
					remove (infile)
				infile = bamfile
			if doIndex:
				cmd = '{{args.picard}} BuildBamIndex %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s"' % (mem, tmpdir, tmpdir, infile, {{idxfile | quote}})
				runcmd (cmd)
			if infile != outfile:
				if path.exists(infile + '.bai'):
					move (infile + '.bai', outfile + '.bai')
				move (infile, outfile)
	
	if not path.exists ({{idxfile | quote}}):
		symlink (outfile, {{idxfile | quote}})
		
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""

"""
@name:
	pBamMarkdup
@description:
	Mark/remove duplicates for bam files
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output bam file
@args:
	`tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)
	`sambamba`         : The path of sambamba. Default: sambamba 
	`picard`           : The path of picard. Default: picard 
	`biobambam_bamsort`: The path of biobambam's bamsort. Default: bamsort 
	`samtools`         : The path of samtools. Default: samtools 
	`bamutil`          : The path of bamutil. Default: bam
	`rmdup`            : Do duplicates removing? Default: False
	- Samtools will anyway remove the duplicates
	`tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
	`nthread`          : Default: 1
	- Not available for samtools and picard
	`params`           : Other parameters for `tool`. Defaut: ""
	`mem`              : The max memory to use. Default: "16G"
	- Unit could be G/g/M/m
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
@requires:
	[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
	[bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)
"""
pBamMarkdup                        = proc (desc = 'Mark/remove duplicates for bam files.')
pBamMarkdup.input                  = "infile:file"
pBamMarkdup.output                 = "outfile:file:{{infile | fn}}.bam"
pBamMarkdup.args.tool              = "biobambam" 
pBamMarkdup.args.sambamba          = "sambamba"
pBamMarkdup.args.picard            = "picard"
pBamMarkdup.args.biobambam_bamsort = "bamsort"
pBamMarkdup.args.samtools          = "samtools"
pBamMarkdup.args.bamutil           = "bam"
pBamMarkdup.args.rmdup             = False
pBamMarkdup.args.tmpdir            = __import__('tempfile').gettempdir() 
pBamMarkdup.args.nthread           = 1
pBamMarkdup.args.params            = ""
pBamMarkdup.args.mem               = "16G"	
pBamMarkdup.args._mem2M            = mem.toM.python
pBamMarkdup.args._memtoJava         = mem.toJava.python
pBamMarkdup.args._runcmd           = runcmd.python
pBamMarkdup.lang                   = "python"
pBamMarkdup.script                 = """
from shutil import move, rmtree
from os import makedirs, path, symlink, remove
from sys import stdout, stderr, exit

{{ args._runcmd }}
{{ args._mem2M }}
{{ args._memtoJava }}

infile    = {{ infile | quote }}
outfile   = {{ outfile | quote }}
tool      = {{ args.tool | quote }}
tmpdir    = {{ args.tmpdir | quote }}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{infile | fn}}.{{#}}")
doRmdup   = {{ args.rmdup }}

if not path.exists (tmpdir):
	makedirs (tmpdir)
try:
	if tool == 'biobambam':
		rmdup = ''
		if doRmdup:
			rmdup = 'rmdup=1 D=/dev/null'
		cmd = '{{args.biobambam_bamsort}} I="%s" O="%s" tmpfile="%s/tmp." %s markthreads={{args.nthread}} {{args.params}}' % (infile, outfile, tmpdir, rmdup)
		runcmd (cmd)
	elif tool == 'sambamba':
		rmdup = ""
		if doRmdup:
			rmdup = "-r"
		cmd = '{{args.sambamba}} markdup %s -t {{args.nthread}} --tmpdir="%s" "%s" "%s"' % (rmdup, tmpdir, infile, outfile)
		runcmd (cmd)
	elif tool == 'samtools':
		cmd = '{{args.samtools}} rmdup "%s" "%s"' % (infile, outfile)
		runcmd (cmd)
	elif tool == 'picard':
		mem = memtoJava({{ args.mem | quote }})
		
		rmdup = ""
		if doRmdup:
			rmdup = "REMOVE_DUPLICATES=true"
		mfile = "/dev/null"
		cmd = '{{args.picard}} MarkDuplicates %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" O="%s" M="%s" %s ' % (mem, tmpdir, tmpdir, infile, outfile, mfile, rmdup)
		runcmd (cmd)
	elif tool == 'bamutil':
		rmdup = ""
		if doRmdup:
			rmdup = "--rmDups"
		cmd = '{{args.bamutil}} dedup --in "%s" --out "%s" %s' % (infile, outfile, rmdup)
		runcmd (cmd)
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""

"""
@name:
	pBamRecal
@description:
	Recalibrate a bam file
@input:
	`infile:file`: The bam file
@brings:
	`infile`: {{infile.orig | bn}}.bai, the index file of bam
@output:
	`outfile:file`: The output bam file
@args:
	`tool`                         : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)
	`gatk`                         : The path of gatk, including java path. Default: `gatk`
	`samtools`                     : The path of samtools. Default: `samtools`
	`bamutil`                      : The path of bamutil. Default: `bam`
	`picard`                       : The path of picard. Default: `picard`
	`paramsRealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
	`paramsIndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
	`paramsBaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
	`paramsPrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
	`params`                       : Other parameters for `bam recab`. Default: ""
	`mem`                          : The max memory to use. Default: "32G"
	`knownSites`                   : The known polymorphic sites to mask out. Default: "" (Required for GATK)
	`ref`                          : The reference file. Required.
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
@requires:
	[gatk](https://software.broadinstitute.org/gatk)
	[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
"""
pBamRecal                                   = proc (desc = 'Recalibrate a bam file.')
pBamRecal.input                             = "infile:file"
pBamRecal.brings                            = {"infile": "{{infile.orig | bn}}.bai"}
pBamRecal.output                            = "outfile:file:{{infile | fn}}.bam, idxfile:file:{{infile | fn}}.bam.bai"
pBamRecal.args.tool                         = "bamutil" 
pBamRecal.args.gatk                         = "gatk" 
pBamRecal.args.samtools                     = "samtools" 
pBamRecal.args.picard                       = "picard" 
pBamRecal.args.bamutil                      = "bam" 
pBamRecal.args.paramsRealignerTargetCreator = ""
pBamRecal.args.paramsIndelRealigner         = ""
pBamRecal.args.paramsBaseRecalibrator       = ""
pBamRecal.args.paramsPrintReads             = ""
pBamRecal.args.params                       = ""
pBamRecal.args.ref                          = ""
pBamRecal.args.tmpdir                       = __import__('tempfile').gettempdir() 
pBamRecal.args.knownSites                   = ""
pBamRecal.args.mem                          = "32G"
pBamRecal.args._memtoJava                   = mem.toJava.python
pBamRecal.args._runcmd                      = runcmd.python
pBamRecal.args._buildArgIndex               = buildArgIndex.python
pBamRecal.beforeCmd                         = checkArgsRef.bash + buildArgsFastaFai.bash + buildArgsFastaDict.bash
pBamRecal.lang                              = "python"
pBamRecal.script                            = """

from os import path, makedirs, remove, symlink
from sys import stderr, exit
from shutil import rmtree, move
from time import sleep

if not path.exists ({{brings.infile | quote}}):
	stderr.write ("Input file '{{infile.orig}}' is not indexed.")
	exit (1)

tool = {{args.tool | quote}}
if not path.exists ({{args.knownSites | quote}}) and tool == 'gatk':
	stderr.write ("knownSites file is required by GATK but is not specified (args.knownSites) or not exists.")
	exit (1)
	
{{ args._runcmd }}
{{ args._memtoJava }}
{{ args._buildArgIndex }}

ref = {{args.ref | quote}}
mem = memtoJava({{ args.mem | quote }})

tmpdir    = {{ args.tmpdir | quote }}
tmpdir    = path.join (tmpdir, "{{proc.id}}.{{infile | fn}}.{{#}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

try:
	if tool == 'gatk':
		intfile = "{{outfile | prefix}}.intervals"
		runcmd ('{{args.gatk}} -T RealignerTargetCreator %s -Djava.io.tmpdir="%s" -R "%s" -I "{{infile}}" -o "%s" {{args.paramsRealignerTargetCreator}}' % (mem, tmpdir, ref, intfile))
		bamfileIr = "{{outfile | prefix}}.ir.bam"
		runcmd ('{{args.gatk}} -T IndelRealigner %s -Djava.io.tmpdir="%s" -R "%s" -I "{{infile}}" -o "%s" -targetIntervals "%s" {{args.paramsIndelRealigner}}' % (mem, tmpdir, ref, bamfileIr, intfile))
		recaltable = "{{outfile | prefix}}.recaltable"
		runcmd ('{{args.gatk}} -T BaseRecalibrator %s -Djava.io.tmpdir="%s" -R "%s" -I "%s" -o "%s" -knownSites "{{args.knownSites}}" {{args.paramsBaseRecalibrator}}' % (mem, tmpdir, ref, bamfileIr, recaltable))
		runcmd ('{{args.gatk}} -T PrintReads %s -Djava.io.tmpdir="%s" -R "%s" -I "%s" -o "{{outfile}}" -BQSR "%s" {{args.paramsPrintReads}}' % (mem, tmpdir, ref, bamfileIr, recaltable))
		remove (bamfileIr)
		move ("{{outfile | prefix}}.bai", "{{idxfile}}")
	elif tool == 'bamutil':
		knownSites = ''
		if "{{args.knownSites}}":
			knownSites = '--dbsnp "{{args.knownSites}}"'
		ref2 = "{{ proc.workdir }}/0/{{ args.ref | bn }}"
		cmd = '{{args.bamutil}} recab --in "{{infile}}" --out "{{outfile}}" --refFile "%s" %s '
		
		refcache = "{{ args.ref | prefix }}-bs.umfa"
		if path.exists (refcache):
			runcmd (cmd % (ref, knownSites))
		else:
			r = buildArgIndex (
				{{#}},
				ref,
				[refcache],
				cmd % (ref, knownSites),
				{{proc.workdir | quote}},
				cmd % (ref2, knownSites)
			)
			if {{#}} != 0:
				runcmd (cmd % (r, knownSites))
		
		cmd = '{{args.samtools}} index "{{outfile}}" "{{idxfile}}"'
		runcmd (cmd)
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)

"""

"""
@name:
	pBamReadGroup
@description:
	Add or replace read groups of a bam file
@input:
	`infile:file`: The bam file
@output:
	`outfile:file`: The output bam file
@args:
	`tool`                         : The tool used. Default: `picard` (picard|bamutil)
	`picard`                       : The path of picard. Default: `picard`
	`bamutil`                      : The path of bamutil. Default: `bam`
	`rg`                           : The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
	- `id` will be parsed from filename with "_LX_" in it if not given
	- `sm` will be parsed from filename
	`params`                       : Other parameters for `tool`. Defaut: ""
	`mem`                          : The max memory to use. Default: "4G"
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
	`tmpdir`                       : The temporary directory. Default: <system tmpdir>
@requires:
	[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
	[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
"""
pBamReadGroup                    = proc (desc = 'Add or replace read groups of a bam file.')
pBamReadGroup.input              = "infile:file"
pBamReadGroup.output             = "outfile:file:{{infile | bn}}"
pBamReadGroup.args.tool          = "bamutil" 
pBamReadGroup.args.picard        = "picard" 
pBamReadGroup.args.bamutil       = "bam" 
pBamReadGroup.args.rg            = {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
pBamReadGroup.args.params        = ""
pBamReadGroup.args.tmpdir        = __import__('tempfile').gettempdir() 
pBamReadGroup.args.mem           = "4G"
pBamReadGroup.args._memtoJava    = mem.toJava.python
pBamReadGroup.args._runcmd       = runcmd.python
pBamReadGroup.lang               = "python"
pBamReadGroup.script             = """
import re
from os import makedirs, path
from shutil import rmtree
from sys import exit, stderr

# detemine default read group
rg = {{ args.rg | json }}
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', "{{outfile | fn}}")
	rg['ID'] = g.group(1) if g else "{{outfile | fn}}.L{{#}}"
if not rg['SM']:
	rg['SM'] = "{{outfile | fn}}"


{{ args._runcmd }}
{{ args._memtoJava }}
mem = memtoJava({{ args.mem | quote }})

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{infile | fn}}.{{#}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

tool = {{args.tool | quote}}
try:
	if tool == 'picard':
		runcmd ('{{args.picard}} AddOrReplaceReadGroups %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="{{infile}}" O="{{outfile}}" RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s' % (mem, tmpdir, tmpdir, rg['ID'], rg['LB'], rg['PL'], rg['PU'], rg['SM']))
	elif tool == 'bamutil':
		runcmd ('{{args.bamutil}} polishBam --RG "@RG\\tID:%s\\t%s" --in "{{infile}}" --out "{{outfile}}"' % (rg['ID'], "\\t".join([k + ":" + v for k,v in rg.items() if k!='ID'])))
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""


"""
@name:
	pBamReorder
@description:
	Reorder a sam/bam file by a given reference file using `picard ReorderSam`
@input:
	`infile:file`: The sam/bam file
@output:
	`outfile:file`: The output bam file
@args:
	`picard`                       : The path of picard. Default: `picard`
	`ref`                          : The reference file. Required
	`params`                       : Other parameters for `picard ReorderSam`. Defaut: ""
	`mem`                          : The max memory to use. Default: "4G"
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
	`tmpdir`                       : The temporary directory. Default: <system tmpdir>
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
"""
pBamReorder                     = proc (desc = 'Reorder a sam/bam file by a given reference.')
pBamReorder.input               = "infile:file"
pBamReorder.output              = "outfile:file:{{infile | bn}}"
pBamReorder.args.picard         = "picard" 
pBamReorder.args.params         = ""
pBamReorder.args.tmpdir         = __import__('tempfile').gettempdir() 
pBamReorder.args.mem            = "4G"
pBamReorder.args.ref            = ""
pBamReorder.args._memtoJava     = mem.toJava.python
pBamReorder.args._runcmd        = runcmd.python
pBamReorder.beforeCmd           = checkArgsRef.bash + buildArgsFastaDict.bash
pBamReorder.lang                = "python"
pBamReorder.script              = """
from os import makedirs, path
from shutil import rmtree
	
tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{infile | fn}}.{{#}}")
if not path.exists (tmpdir): makedirs (tmpdir)
	
{{ args._runcmd }}
{{ args._memtoJava }}
mem = memtoJava({{ args.mem | quote }})
ref = {{args.ref | quote}}	
try:
	runcmd ('{{args.picard}} ReorderSam %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="{{infile}}" O="{{outfile}}" R="%s" {{args.params}}' % (mem, tmpdir, tmpdir, ref))
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""


"""
@name:
	pBamMerge
@description:
	Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.
@input:
	`inlist:file`: The directory containing sam/bam files to be merged
@output:
	`outfile:file`: The merged bam file
@args:
	`tool`     : The tool used to merge. Default: bamutil (picard|samtools|sambamba)
	`picard`   : The path of picard. Default: `picard`
	`bamutil`  : The path of bamutil. Default: `bam`
	`samtools` : The path of samtools. Default: `samtools`
	`sambamba` : The path of sambamba. Default: `sambamba`
	`params`   : Other parameters for `tool`. Defaut: ""
	`mem`      : The max memory to use. Default: "4G"
	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
	`tmpdir`   : The temporary directory. Default: <system tmpdir>
	`nthread`  : # threads to use. Default: 1
	- For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
"""
pBamMerge                     = proc (desc = 'Merges multiple SAM and/or BAM sorted files into a single file.')
pBamMerge.input               = "inlist:file"
pBamMerge.output              = "outfile:file:{{inlist | fn | lambda x: x + '_merged'}}.bam"
pBamMerge.args.tool           = "picard" 
pBamMerge.args.picard         = "picard" 
pBamMerge.args.bamutil        = "bam" 
pBamMerge.args.samtools       = "samtools" 
pBamMerge.args.sambamba       = "sambamba" 
pBamMerge.args.params         = ""
pBamMerge.args.tmpdir         = __import__('tempfile').gettempdir() 
pBamMerge.args.nthread        = 1
pBamMerge.args.mem            = "4G"
pBamMerge.args._memtoJava     = mem.toJava.python
pBamMerge.args._runcmd        = runcmd.python
pBamMerge.lang                = "python"
pBamMerge.script              = """

from os import makedirs, path
from shutil import rmtree

tmpdir    = path.join ("{{ args.tmpdir }}", "{{proc.id}}.{{inlist | fn}}.{{#}}")
if not path.exists (tmpdir): makedirs (tmpdir)
	
{{ args._runcmd }}
{{ args._memtoJava }}
mem = memtoJava({{ args.mem | quote }})

tool = {{args.tool | quote}}
try:
	if tool == 'picard':
		infiles = {{ inlist | readlines | json }}
		infiles = " ".join(['I="%s"' % x for x in infiles])
		thr = "USE_THREADING=true" if {{args.nthread}} > 1 else "USE_THREADING=false"
		cmd = '{{args.picard}} MergeSamFiles %s -Djava.io.tmpdir="%s" TMP_DIR="%s" %s O="{{outfile}}" AS=true %s' % (mem, tmpdir, tmpdir, infiles, thr)
		runcmd (cmd)
	elif tool == 'bamutil':
		infiles = {{ inlist | readlines | json }}
		infiles = " ".join(['-i "%s"' % x for x in infiles])
		cmd = '{{args.bamutil}} mergeBam %s -o "{{outfile}}"' % (infiles)
		runcmd (cmd)
	elif tool == 'samtools':
		cmd = '{{args.samtools}} merge -@ {{args.nthread}} -O bam -b "{{inlist}}" "{{outfile}}"'
		runcmd (cmd)
	elif tool == 'sambamba':
		infiles = {{ inlist | readlines | json }}
		infiles = ' '.join(['"'+infile+'"' for infile in infiles])
		cmd = '{{args.sambamba}} merge -t {{args.nthread}} "{{outfile}}" %s' % (infiles)
		runcmd (cmd)
except Exception as ex:		
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""


"""
@name:
	pBam2Gmut
@description:
	Call germline (snps and indels) from a call-ready bam file.
@input:
	`infile:file`: The input bam file
@brings:
	`infile`: `{{infile | bn}}.bai`, the bam index file
@output:
	`outfile:file`: The vcf file containing the mutations
@args:
	`tool`:         The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)
	`gatk`:         The path of gatk. Default: gatk
	`vardict`:      The path of vardict. Default: vardict
	`snvsniffer`:   The path of snvsniffer. Default: SNVSniffer
	`samtools`:     The path of samtools. Default: samtools (used to generate reference index)
	`platypus`:     The path of platypus. Default: platypus
	`strelka`:      The path of strelka. Default: configureStrelkaGermlineWorkflow.py
	`configParams`: The params for `strelka` configuration. Default: ""
	`picard`:       The path of picard. Default: picard
	`mem`:          The memory to be used. Default: 32G
	- will be converted to -Xms4G -Xmx32G for java programs
	`ref`:          The reference file. Required.
	`gz`:           Gzip output file? Default: False
	`tmpdir`:       The temporary directory. Default: <system tmpdir>
	`params`:       Other params for `tool`. Default: ""
@requires:
	[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
	[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
	[vardict](https://github.com/AstraZeneca-NGS/VarDict)
	[snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
	[platypus](http://www.well.ox.ac.uk/platypus)
	[strelka@2.7.1+](https://github.com/Illumina/strelka)
"""
pBam2Gmut                     = proc (desc = 'Call germline (snps and indels) from a call-ready bam file.')
pBam2Gmut.input               = "infile:file"
pBam2Gmut.brings              = {"infile": "{{infile.orig | bn}}.bai"}
pBam2Gmut.output              = "outfile:file:{{infile | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pBam2Gmut.lang                = "python"
pBam2Gmut.args.tool           = "gatk"
pBam2Gmut.args.gatk           = "gatk"
pBam2Gmut.args.vardict        = "vardict"
pBam2Gmut.args.snvsniffer     = "SNVSniffer"
pBam2Gmut.args.samtools       = "samtools" # required by SNVSniffer to generate a bam header file
pBam2Gmut.args.platypus       = "platypus"
pBam2Gmut.args.picard         = "picard"
pBam2Gmut.args.strelka        = "configureStrelkaGermlineWorkflow.py"
pBam2Gmut.args.mem            = "24G"
pBam2Gmut.args.ref            = ""
pBam2Gmut.args.tmpdir         = __import__('tempfile').gettempdir() 
pBam2Gmut.args.configParams   = '' # only for strelka
pBam2Gmut.args.params         = ""
pBam2Gmut.args.gz             = False
pBam2Gmut.args.nthread        = 1 # for gatk and platypus
pBam2Gmut.args._memtoJava     = mem.toJava.python
pBam2Gmut.args._runcmd        = runcmd.python
pBam2Gmut.beforeCmd           = checkArgsRef.bash + buildArgsFastaFai.bash + buildArgsFastaDict.bash
pBam2Gmut.script              = """
from os import path, makedirs
from shutil import rmtree
from sys import stderr, exit

if not path.exists ("{{brings.infile}}"):
	stderr.write ("Input file '{{infile.orig}}' is not indexed.")
	exit (1)
	
{{ args._runcmd }}
{{ args._memtoJava }}

tmpdir    = path.join ("{{ args.tmpdir}}", "{{proc.id}}.{{infile | fn}}.{{#}}")
if not path.exists (tmpdir): makedirs (tmpdir)

ref     = {{args.ref | quote}}
tool    = {{args.tool | quote}}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{outfile | quote}}
if gz:	outfile = outfile[:-3]
try:
	if tool == 'gatk':
		mem = memtoJava ({{args.mem | quote}})
		cmd = '{{args.gatk}} -T HaplotypeCaller %s -Djava.io.tmpdir="%s" -R "%s" -I "{{infile}}" -o "%s" -nct {{args.nthread}} {{args.params}}' % (mem, tmpdir, ref, outfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'vardict':	
		cmd = '{{args.vardict}} -v -G "%s" -b "{{infile}}" {{args.params}} > "%s"' % (ref, outfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'snvsniffer':
		hfile = "{{job.outdir}}/{{infile | bn}}.header"
		cmd = '{{args.samtools}} view -H "{{infile}}" > "%s"' % hfile
		runcmd (cmd)
		cmd = '{{args.snvsniffer}} snp -g "%s" -o "%s" {{args.params}} "%s" "{{infile}}"' % (ref, outfile, hfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'platypus':
		cmd = '{{args.platypus}} callVariants --refFile="%s" --bamFiles="{{infile}}" --nCPU={{args.nthread}} --output="%s" --logFileName="%s"' % (ref, outfile, outfile + '.log')
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'strelka':
		# config
		cmd = '{{args.strelka}} --bam="{{infile}}" --referenceFasta="%s" --runDir="{{job.outdir}}" {{args.configParams}}' % (ref)
		runcmd (cmd)
		# run
		cmd = '{{job.outdir}}/runWorkflow.py -m local -j "{{args.nthread}}" -g {{args.mem | lambda x: int(x[:-1]) if x.lower().endswith("g") else int(int(x[:-1])/1024) }} {{args.params}}'
		runcmd (cmd)
		ofile = "{{job.outdir}}/results/variants/genome.S1.vcf.gz"
		runcmd ('mv "%s" "%s.gz"' % (ofile, outfile))
		if not gz: runcmd ('gunzip "%s.gz"' % outfile)
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""

pBamPair2Smut                     = proc (desc = 'Call somatic mutations from tumor-normal bam pair.')
pBamPair2Smut.input               = "tumor:file, normal:file"
pBamPair2Smut.brings              = {"tumor": "{{tumor.orig | bn}}.bai", "normal": "{{normal.orig | bn}}.bai"}
pBamPair2Smut.output              = "outfile:file:{{tumor | fn | fn}}-{{normal | fn | fn}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}"
pBamPair2Smut.args.tool           = 'gatk'
pBamPair2Smut.args.gatk           = 'gatk' # required for strelka
pBamPair2Smut.args.somaticsniper  = 'bam-somaticsniper'
pBamPair2Smut.args.strelka        = 'configureStrelkaSomaticWorkflow.py' # @2.7.1
pBamPair2Smut.args.snvsniffer     = 'SNVSniffer'
pBamPair2Smut.args.virmid         = 'virmid'
pBamPair2Smut.args.vardict        = 'vardict'
pBamPair2Smut.args.samtools       = 'samtools'
pBamPair2Smut.args.picard         = 'picard'
pBamPair2Smut.args.configParams   = '' # only for strelka
pBamPair2Smut.args.params         = '' 
pBamPair2Smut.args.mem            = '24G'
pBamPair2Smut.args.ref            = ''
pBamPair2Smut.args.gz             = False
pBamPair2Smut.args.nthread        = 1
pBamPair2Smut.args.tmpdir         = __import__('tempfile').gettempdir() 
pBamPair2Smut.args._memtoJava     = mem.toJava.python
pBamPair2Smut.args._runcmd        = runcmd.python
pBamPair2Smut.beforeCmd           = checkArgsRef.bash + buildArgsFastaFai.bash + buildArgsFastaDict.bash
pBamPair2Smut.lang                = 'python'
pBamPair2Smut.script              = """
from os import path, makedirs
from shutil import rmtree
from sys import stderr, exit

if not path.exists ("{{brings.tumor}}"):
	stderr.write ("Input file '{{tumor.orig}}' is not indexed.")
	exit (1)
if not path.exists ("{{brings.normal}}"):
	stderr.write ("Input file '{{normal.orig}}' is not indexed.")
	exit (1)

if not path.exists ({{args.ref | quote}}):
	stderr.write ("Reference file not specified")
	exit (1)

{{ args._runcmd }}
{{ args._memtoJava }}

tmpdir    = path.join ("{{args.tmpdir}}", "{{proc.id}}.{{tumor | fn | fn}}.{{normal | fn | fn}}.{{#}}")
if not path.exists (tmpdir):
	makedirs (tmpdir)

ref     = {{args.ref | quote}}
tool    = {{args.tool | quote}}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{outfile | quote}}
if gz:	outfile = outfile[:-3]
try:
	if tool == 'gatk':
		intvfile = "{{job.outdir}}/interval.list"
		cmd = '{{args.samtools}} idxstats "{{tumor}}" | head -n -1 | cut -f1 > "%s"' % (intvfile)
		runcmd (cmd)
		mem = memtoJava ({{args.mem | quote}})
		cmd = '{{args.gatk}} -T MuTect2 %s -Djava.io.tmpdir="%s" -R "%s" -I:tumor "{{tumor}}" -I:normal "{{normal}}" -o "%s" -nct {{args.nthread}} -L "%s" {{args.params}}' % (mem, tmpdir, ref, outfile, intvfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'somaticsniper':	
		cmd = '{{args.somaticsniper}} -f "%s" -F vcf "{{tumor}}" "{{normal}}" "%s"' % (ref, outfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'snvsniffer':
		theader = "{{job.outdir}}/{{tumor | bn}}.header"
		cmd = '{{args.samtools}} view -H "{{tumor}}" > "%s"' % theader
		runcmd (cmd)
		nheader = "{{job.outdir}}/{{normal | bn}}.header"
		cmd = '{{args.samtools}} view -H "{{normal}}" > "%s"' % nheader
		runcmd (cmd)
		cmd = '{{args.snvsniffer}} somatic -g "%s" -o "%s" {{args.params}} "%s" "%s" "{{tumor}}" "{{normal}}"' % (ref, outfile, theader, nheader)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'strelka':
		# config
		cmd = '{{args.strelka}} --normalBam="{{normal}}" --tumorBam="{{tumor}}" --referenceFasta="%s" --runDir="{{job.outdir}}" {{args.configParams}}' % (ref)
		runcmd (cmd)
		# run
		cmd = '{{job.outdir}}/runWorkflow.py -m local -j {{args.nthread}} -g {{args.mem | lambda x: int(x[:-1]) if x.lower().endswith("g") else int(int(x[:-1])/1024) }} {{args.params}}'
		runcmd (cmd)
		# merge
		mem = memtoJava ({{args.mem | quote}})
		cmd = '{{args.gatk}} -T CombineVariants %s -Djava.io.tmpdir="%s" -R "%s" -V:SNV "{{job.outdir}}/results/variants/somatic.snvs.vcf.gz" -V:INDEL "{{job.outdir}}/results/variants/somatic.indels.vcf.gz" -o "%s" --genotypemergeoption UNIQUIFY' % (mem, tmpdir, ref, outfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'virmid':
		mem = memtoJava ({{args.mem | quote}})
		cmd = '{{args.virmid}} %s -Djava.io.tmpdir="%s" -R "%s" -D "{{tumor}}" -N "{{normal}}" -w "{{job.outdir}}" {{args.params}}' % (mem, tmpdir, ref)
		runcmd (cmd)
		runcmd ('mv "{{job.outdir}}/*.virmid.som.passed.vcf" "%s"' % outfile)
		if gz:	runcmd ('gzip "%s"' % (outfile))
	elif tool == 'vardict':
		cmd = '{{args.vardict}} -v -G "%s" -b "{{tumor}}|{{normal}}" {{args.params}} > "%s"' % (ref, outfile)
		runcmd (cmd)
		if gz:	runcmd ('gzip "%s"' % (outfile))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
"""

"""
@name:
	pBam2Cnv
@description:
	Detect copy number variation from bam files.
@input:
	`input:file`: The bam file
@brings:
	`infile`: "{{infile | bn}}.bai" The bam index file
@output:
	`outfile:file`: The output vcf file
	`outdir`: The output directory containing other result files
@args:
	`gz`                    : Whether to gzip the output vcf file. Default: False
	`tool`                  : The tool used to call cnv. Default: 'cnvkit'
	`cnvnator`              : The path of cnvnator. Default: 'cnvnator'
	`cnvnator2vcf`          : The path of cnvnator2VCF. Default: 'cnvnator2VCF.pl'
	`cnvkit`                : The path of cnvkit. Default: 'cnvkit.py'
	`wandy`                 : Tha path of Wandy. Default: 'Wandy'. A `tool.info` file should be with the executable file.
	`ref`                   : The reference file. Required by cnvkit to generate access file. Default: ''
	`cnvkitAccessParams`    : The params for cnvkit access command. Default: '-s 5000'
	`cnvkitTargetParams`    : The params for cnvkit target command. Default: '--split --short-names'
	`cnvkitCoverageParams`  : The params for cnvkit coverage command. Default: ''
	`cnvkitReferenceParams` : The params for cnvkit reference command. Default: '--no-edge'
	`cnvkitFixParams`       : The params for cnvkit fix command. Default: '--no-edge'
	`cnvkitSegmentParams`   : The params for cnvkit segment command. Default: ''
	`cnvkitCallParams`      : The params for cnvkit call command. Default: ''
	`cnvkitPlotParams`      : The params for cnvkit plot command. Default: ''
	`cnvkitBreaksParams`    : The params for cnvkit breaks command. Default: ''
	`cnvkitGainlossParams`  : The params for cnvkit gainloss command. Default: ''
	`cnvkitMetricsParams`   : The params for cnvkit metrics command. Default: ''
	`cnvkitSegmetricsParams`: The params for cnvkit segmetrics command. Default: '--iqr'
	`cnvkitExportParams`    : The params for cnvkit export command. Default: ''
	`cnvkitScatterParams`   : The params for cnvkit scatter command. Default: [''] # multiple scatter plots
	`cnvkitHeatmapParams`   : The params for cnvkit heatmap command. Default: [''] # multiple heatmap plots
	`cnvkitDiagramParams`   : The params for cnvkit diagram command. Default: ''
	`cnvkitReport`          : Generate cnvkit reports? Default: True
	`cnvkitPlot`            : Generate cnvkit plots? Default: True
	`cnvnatorBinsize`       : Bin size for cnvnator. Default: 100
	`cnvnatorGenome`        : Genome for cnvnator. Default: 'hg19'. (NCBI36, hg18, GRCh37, hg19)
	`params`                : The params for `tool`. Default: '-t 1' # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa
	`mem`                   : The memory used. Default: '20G' # only for wandy
	`nthread`               : The # threads to use. Default: 1	 # only for cnvkit
@requires:
	[`cnvkit`](http://cnvkit.readthedocs.io/en/stable/index.html)
	[`cnvnator`](https://github.com/abyzovlab/CNVnator)
	`wandy`: Inside cnv caller
"""
pBam2Cnv = proc (desc = 'Detect copy number variation from bam files.')
pBam2Cnv.input = 'infile:file'
pBam2Cnv.brings = {"infile": "{{infile | bn}}.bai"}
pBam2Cnv.output = "outfile:file:{{infile | fn}}.{{args.tool}}/{{infile | fn}}.{{args.tool}}.vcf{{args.gz | lambda x: '.gz' if x else ''}}, outdir:dir:{{infile | fn}}.{{args.tool}}"
pBam2Cnv.args.gz = False
pBam2Cnv.args.tool = 'cnvkit'
pBam2Cnv.args.cnvnator = 'cnvnator'
pBam2Cnv.args.cnvnator2vcf = 'cnvnator2VCF.pl'
pBam2Cnv.args.cnvkit = 'cnvkit.py'
pBam2Cnv.args.wandy = 'Wandy'
pBam2Cnv.args.ref = ''
pBam2Cnv.args.cnvkitAccessParams = '-s 5000'
pBam2Cnv.args.cnvkitTargetParams = '--split --short-names'
pBam2Cnv.args.cnvkitCoverageParams = ''
pBam2Cnv.args.cnvkitReferenceParams = '--no-edge'
pBam2Cnv.args.cnvkitFixParams = '--no-edge'
pBam2Cnv.args.cnvkitSegmentParams = ''
pBam2Cnv.args.cnvkitCallParams = ''
pBam2Cnv.args.cnvkitPlotParams = ''
pBam2Cnv.args.cnvkitBreaksParams = ''
pBam2Cnv.args.cnvkitGainlossParams = ''
pBam2Cnv.args.cnvkitMetricsParams = ''
pBam2Cnv.args.cnvkitSegmetricsParams = '--iqr'
pBam2Cnv.args.cnvkitExportParams = ''
pBam2Cnv.args.cnvkitScatterParams = [''] # multiple scatter plots
pBam2Cnv.args.cnvkitHeatmapParams = [''] # multiple heatmap plots
pBam2Cnv.args.cnvkitDiagramParams = ''
pBam2Cnv.args.cnvkitReport = True
pBam2Cnv.args.cnvkitPlot = True
pBam2Cnv.args.cnvnatorBinsize = 100
pBam2Cnv.args.cnvnatorGenome = 'hg19'
pBam2Cnv.args.params = '-t 1' # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa
pBam2Cnv.args.mem = '20G' # for wandy
pBam2Cnv.args.nthread = 1 # for cnvkit
pBam2Cnv.args._polling0 = polling0.python
pBam2Cnv.args._pollingAll = pollingAll.python
pBam2Cnv.args._runcmd = runcmd.python
pBam2Cnv.beforeCmd = """
if [[ ! -e "{{args.ref}}" && "{{args.tool}}" == "cnvkit" ]]; then
	echo "No reference file specified." 1>&2
	exit 1
fi
if [[ "{{args.tool}}" == "cnvkit" && {{proc.forks}} -lt {{proc.length}} ]]; then
	echo "Cnvkit requires all jobs run simultaneously (proc.forks >= # jobs)." 1>&2
	echo "Because it needs all coverage files to generate reference coverage file." 1>&2
	exit 1
fi
"""
pBam2Cnv.lang   = "python"
pBam2Cnv.script = """

from os import path, remove, makedirs
from sys import stderr, exit
from time import sleep
from datetime import date

if not path.exists ("{{brings.infile}}"):
	stderr.write ("Input file '{{infile.orig}}' is not indexed.")
	exit (1)

{{ args._runcmd }}
{{ args._polling0 }}
{{ args._pollingAll }}
gz      = {{args.gz | lambda x: bool(x)}}
outfile = {{outfile | quote}}
if gz:	outfile = outfile[:-3]
	
tool = {{args.tool | quote}}
if tool == 'cnvkit':
	targetDone    = "{{proc.workdir}}/0/output/cnvkit_target.done"
	referenceDone = "{{proc.workdir}}/0/output/cnvkit_reference.done"
	targetCov     = "{{outdir}}/{{infile | fn}}.targetcov.cnn"
	accessfile    = '{{proc.workdir}}/0/output/cnvkit_access.bed'
	targetfile    = '{{proc.workdir}}/0/output/cnvkit_targets.bed'
	refcnn        = '{{proc.workdir}}/0/output/reference.cnn'
	fixedCnr      = "{{outdir}}/{{infile | fn}}.cnr"
	segfile       = "{{outdir}}/{{infile | fn}}.cns"
	callfile      = "{{outdir}}/{{infile | fn}}.call.cns"
	# report files
	breaksfile    = "{{outdir}}/{{infile | fn}}.breaks.txt"
	gainlossfile  = "{{outdir}}/{{infile | fn}}.gainloss.txt"
	metricsfile   = "{{outdir}}/{{infile | fn}}.metrics.txt"
	segmetricsfile= "{{outdir}}/{{infile | fn}}.segmetrics.txt"
	openblas_nthr = "export OPENBLAS_NUM_THREADS={{args.nthread}}; export OMP_NUM_THREADS={{args.nthread}}; export NUMEXPR_NUM_THREADS={{args.nthread}}; export MKL_NUM_THREADS={{args.nthread}}; "
	
	cmd1 = '%s {{args.cnvkit}} access "{{args.ref}}" -o "%s" {{args.cnvkitAccessParams}}' % (openblas_nthr, accessfile)
	cmd2 = '{{args.cnvkit}} target "%s" -o "%s" {{args.cnvkitTargetParams}}' % (accessfile, targetfile)
	polling0 ({{#}}, cmd1 + '; ' + cmd2, targetDone, t = 60)
	
	cmd = '%s {{args.cnvkit}} coverage "{{infile}}" "%s" -p {{args.nthread}} -o "%s" {{args.cnvkitCoverageParams}}' % (openblas_nthr, targetfile, targetCov)
	pollingAll ({{proc.workdir | quote}}, {{proc.length}}, {{#}}, cmd, "cnvkit_coverage.done")
	
	cmd = '%s {{args.cnvkit}} reference {{proc.workdir}}/*/output/*/*.targetcov.cnn {{args.cnvkitReferenceParams}} -o "%s" -f "{{args.ref}}"' % (openblas_nthr, refcnn)
	polling0 ({{#}}, cmd, referenceDone)
	
	mtfile = "{{outdir}}/cnvkit_mt"
	open(mtfile, 'w').close()
	cmd = '%s {{args.cnvkit}} fix "%s" "%s" "%s" {{args.cnvkitFixParams}} -o "%s"' % (openblas_nthr, targetCov, mtfile, refcnn, fixedCnr)
	runcmd (cmd)
	if path.getsize(fixedCnr) < 60:
		open(segfile, 'w').write('chromosome	start	end	gene	log2	depth	probes	weight\\n')
	else:
		cmd = '%s {{args.cnvkit}} segment -o "%s" -p {{args.nthread}} {{args.cnvkitSegmentParams}} "%s"' % (openblas_nthr, segfile, fixedCnr)
		runcmd (cmd)
	if path.getsize(segfile) < 60:
		open(callfile, 'w').write('chromosome	start	end	gene	log2	cn	depth	probes	weight\\n')
	else:
		cmd = '%s {{args.cnvkit}} call -o "%s" {{args.cnvkitCallParams}} "%s"' % (openblas_nthr, callfile, segfile)
		runcmd (cmd)
	
	# reports
	if {{args.cnvkitReport | lambda x: bool(x)}}:
		cmd = '%s {{args.cnvkit}} breaks "%s" "%s" -o "%s" {{args.cnvkitBreaksParams}}' % (openblas_nthr, fixedCnr, callfile, breaksfile)
		runcmd (cmd)
		cmd = '%s {{args.cnvkit}} gainloss "%s" -s "%s" -o "%s" {{args.cnvkitGainlossParams}}' % (openblas_nthr, fixedCnr, callfile, gainlossfile)
		runcmd (cmd)
		cmd = '%s {{args.cnvkit}} metrics "%s" -s "%s" -o "%s" {{args.cnvkitMetricsParams}}' % (openblas_nthr, fixedCnr, callfile, metricsfile)
		runcmd (cmd)
		cmd = '%s {{args.cnvkit}} segmetrics "%s" -s "%s" -o "%s" {{args.cnvkitSegmetricsParams}}' % (openblas_nthr, fixedCnr, callfile, segmetricsfile)
		runcmd (cmd)
	
	if {{args.cnvkitPlot | lambda x: bool(x)}}:
	# plots
		#scatter plots
		for i, param in enumerate({{args.cnvkitScatterParams | json}}):
			cmd = '%s {{args.cnvkit}} scatter "%s" -s "%s" %s -o "{{outdir}}/{{infile | fn}}.scatter%s.pdf"' % (openblas_nthr, fixedCnr, callfile, param, str(i))
			runcmd (cmd)
		cmd = '%s {{args.cnvkit}} diagram "%s" -s "%s" {{args.cnvkitDiagramParams}} -o {{outdir}}/{{infile | fn}}.diagram.pdf' % (openblas_nthr, fixedCnr, callfile)
		runcmd (cmd)
		
		pollingAll ({{proc.workdir | quote}}, {{proc.length}}, {{#}}, 'if [[ ! -e "%s" ]]; then sleep 1; fi' % callfile, "cnvkit_call.done")
		if {{#}} == 0:
			for i, param in enumerate({{args.cnvkitHeatmapParams | json}}):
				cmd = '%s {{args.cnvkit}} heatmap "{{proc.workdir}}"/*/output/*/*.call.cns %s -o "{{outdir}}/{{infile | fn}}.heatmap%s.pdf"' % (openblas_nthr, param, str(i))
				runcmd (cmd)
	
	# to vcf
	if path.getsize(callfile) < 60: # no data generated
		with open (outfile, 'w') as fout:
			fout.write('##fileformat=VCFv4.2\\n')
			fout.write('##fileDate=%s\\n' % str(date.today()).replace('-', ''))
			fout.write('##source=CNVkit\\n')
			fout.write('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">\\n')
			fout.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\\n')
			fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\\n')
			fout.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\\n')
			fout.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\\n')
			fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\\n')
			fout.write('##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description="Fold change">\\n')
			fout.write('##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">\\n')
			fout.write('##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">\\n')
			fout.write('##ALT=<ID=DEL,Description="Deletion">\\n')
			fout.write('##ALT=<ID=DUP,Description="Duplication">\\n')
			fout.write('##ALT=<ID=CNV,Description="Copy number variable region">\\n')
			fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\\n')
			fout.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">\\n')
			fout.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\\n')
			fout.write('##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">\\n')
			fout.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{{infile | fn}}\\n')
	else:
		cmd = '%s {{args.cnvkit}} export vcf "%s" -o "%s" {{args.cnvkitExportParams}}' % (openblas_nthr, callfile, outfile)
		runcmd (cmd)
	
elif tool == 'cnvnator':
	rootfile = '{{outdir}}/{{infile | fn}}.root'
	callfile = '{{outdir}}/{{infile | fn}}.cnvnator'
	cmd = '{{args.cnvnator}} -root "%s" -genome {{args.cnvnatorGenome}} -tree "{{infile}}" ' % (rootfile)
	runcmd (cmd)
	cmd = '{{args.cnvnator}} -root "%s" -genome {{args.cnvnatorGenome}} -his {{args.cnvnatorBinsize}}' % (rootfile)
	runcmd (cmd)
	cmd = '{{args.cnvnator}} -root "%s" -stat {{args.cnvnatorBinsize}}' % (rootfile)
	runcmd (cmd)
	cmd = '{{args.cnvnator}} -root "%s" -partition {{args.cnvnatorBinsize}}' % (rootfile)
	runcmd (cmd)
	cmd = '{{args.cnvnator}} -root "%s" -call {{args.cnvnatorBinsize}} > "%s"' % (rootfile, callfile)
	runcmd (cmd)
	cmd = 'cd "{{outdir}}"; {{args.cnvnator2vcf}} "%s" > "%s"' % (path.basename(callfile), outfile)
	runcmd (cmd)

elif tool == 'wandy':
	# get Wandy tool.info
	from distutils.spawn import find_executable
	toolinfo   = path.join (path.dirname(find_executable("{{args.wandy}}")), "tool.info")
	retdir     = "{{outdir}}"
	#if not path.exists(retdir):
	#	makedirs(retdir)
	myToolinfo = path.join (retdir, "tool.info")
	if not path.exists(toolinfo):
		stderr.write('Cannot find tool.info in wandy source directory.')
		exit(1)
	with open (toolinfo) as fin, open(myToolinfo, 'w') as fout:
		toolinfos  = [line.strip() for line in fin if not line.strip().startswith('WANDY_MEM') and not line.strip().startswith('WANDY_PARALLELED_MODE')]
		toolinfos.append ('WANDY_MEM={{args.mem}}')
		toolinfos.append ('WANDY_PARALLELED_MODE=0')
		fout.write ("\\n".join(toolinfos) + "\\n")
		
	cmd = 'cd "%s"; {{args.wandy}} -i "{{infile}}" {{args.params}}' % (retdir)
	runcmd (cmd)
	
	# TODO: convert to vcf
	open({{outfile | quote}}, 'w').close()
	
if gz:	runcmd ('gzip "%s"' % (outfile))
"""

pBamStats = proc (desc = 'Get read depth from bam files.')
pBamStats.input = 'infile:file'
pBamStats.output = 'outfile:file:{{infile | fn}}/{{infile | fn}}.stat.txt, outdir:dir:{{infile | fn}}'
pBamStats.args.tool = 'bamstats'
pBamStats.args.bamstats = 'bamstats'
pBamStats.args.params = ''
pBamStats.args.mem = '16G'
pBamStats.args.plot = True
pBamStats.args._runcmd = runcmd.r
pBamStats.args._memtoJava = mem.toJava.r
pBamStats.args._pollingAll = pollingAll.r
pBamStats.beforeCmd = """
if [[ "{{args.plot | Rbool}}" == "TRUE" && {{proc.forks}} -lt {{proc.length}} ]]; then
	echo "Plots can only be done with all jobs run simultaneously (proc.forks >= # jobs)." 1>&2
	exit 1
fi
"""
pBamStats.lang = 'Rscript'
pBamStats.script = """
{{args._memtoJava}}

cmd = '{{args.bamstats}} -i "{{infile}}" -o "{{outfile}}" {{args.params}}'

if ({{args.plot | Rbool}}) {
	{{args._pollingAll}}
	pollingAll ({{proc.workdir | quote}}, {{proc.length}}, {{#}}, cmd, "bamstats.done")
	
	##### start plotting
	if ({{#}} == 0) {
	
		bsfiles = Sys.glob("{{proc.workdir}}/*/output/*/*.stat.txt")
		means   = matrix(ncol=1, nrow=length(bsfiles))
		chrs    = NULL
		
		rnames  = make.unique(unlist(lapply(bsfiles, function(x){ x=basename(x); return (substr(x, 1, nchar(x)-9)) })), sep='_')
		rownames(means) = rnames
		colnames(means) = c("Average coverage")
		for (i in 1:length(bsfiles)) {
			bsfile = bsfiles[i]
			write (paste("Reading", bsfile, "...", sep=" "), stderr())
			sample = rnames[i]
			stat   = read.table (bsfile, sep="", header=T, check.names=F, row.names=1)
			#stat   = stat[chrs2, ]
			stat[, "N"] = as.numeric(gsub(",", "", stat[, "N"]))
			means[sample, 1] = sum(stat[, "N"] * stat[, "mean"])/sum(stat[, "N"])
			col2in = stat[, "mean", drop=F]
			colnames(col2in) = sample
			if (is.null(chrs)) {
				chrs = col2in
			} else {
				chrs = cbind(chrs, col2in)
			}
		}
		# plot average coverage
		plotFreq = function (obj, figure, xlab, ylab="Frequency") {
			png (file=figure)
			h = hist (obj, freq=T, xlab=xlab, ylab=ylab, col="gray", main=paste(xlab, "distribution", sep=" "), axes=F)
			minb = min(h$breaks)
			maxb = max(h$breaks)
			maxc = max(h$counts)
			lenb = length(h$breaks)
			stpb = (maxb-minb)/(lenb-1)
			axis(1, pos=0, labels=T, at=seq(minb,maxb,stpb))
			lab0 = floor(log10(maxc))
			stpc = ceiling(maxc/(10**lab0)) * (10 ** (lab0-1))
			axis(2, pos=minb, labels=T, at=seq(0, maxc, stpc))
			dev.off()
		}
		write ("Plotting average coverages ...", stderr())
		write.table (means, "{{outdir}}/avgCoverage.txt", quote=F, sep="\\t")
		plotFreq (means, "{{outdir}}/avgCoverage.png", xlab="Average coverage")

		# plot chromosomes
		write ("Plotting chromosome coverages ...", stderr())
		png ("{{outdir}}/chrCoverage.png")
		#colnames(chrs) = NULL # just show index
		boxplot(t(chrs), ylab="Coverage", las=2)
		dev.off()
	}
	
} else {
	runcmd (cmd)
}
"""

"""
@name:
	pBam2FastqPE
@description:
	Convert sam/bam files to pair-end fastq files.
@input:
	`infile:file`: The sam/bam file. 
		- Sam files only available for biobambam, picard
@output:
	`fqfile1:file`: The 1st match of paired reads
	`fqfile2:file`: The 2nd match of paired reads
@args:
	`tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
	`biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
	`bedtools`            : The path of bedtools. Default: bedtools
	`samtools`            : The path of samtools. Default: samtools
	`picard`              : The path of picard. Default: picard
	`mem`                 : The memory to be used by picard. Default: 8G
	`gz`                  : Whether gzip the output files. Default: True
	`params`:             : Other params for `tool`. Default: ''
	`tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
	[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
"""
pBam2FastqPE                           = proc (desc = 'Convert bam files to pair-end fastq files.')
pBam2FastqPE.input                     = "infile:file"
pBam2FastqPE.output                    = [
	"fqfile1:file:{{ infile | fn | fn | fn }}_1.fq{{args.gz | lambda x: '.gz' if x else ''}}", 
	"fqfile2:file:{{ infile | fn | fn | fn }}_2.fq{{args.gz | lambda x: '.gz' if x else ''}}"
]
pBam2FastqPE.args.tool                 = 'bedtools'
pBam2FastqPE.args.biobambam_bamtofastq = 'bamtofastq'
pBam2FastqPE.args.bedtools             = 'bedtools'
pBam2FastqPE.args.samtools             = 'samtools'
pBam2FastqPE.args.picard               = 'picard'
pBam2FastqPE.args.mem                  = '8G' # only for picard
pBam2FastqPE.args.gz                   = True
pBam2FastqPE.args.params               = ''
pBam2FastqPE.args.tmpdir               = __import__('tempfile').gettempdir() 
pBam2FastqPE.args._runcmd              = runcmd.python
pBam2FastqPE.args._memtoJava           = mem.toJava.python
pBam2FastqPE.lang                      = 'python'
pBam2FastqPE.script                    = """
from os import makedirs, path

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

infile  = {{infile | quote}}
fqfile1 = {{fqfile1 | quote}}
fqfile2 = {{fqfile2 | quote}}
if {{args.gz | bool}}:
	fqfile1 = fqfile1[:-3]
	fqfile2 = fqfile2[:-3]

{{args._runcmd}}
{{args._memtoJava}}

params  = {{args.params | quote}}
tool    = {{args.tool | quote}}
try:
	if tool == 'biobambam':
		params += ' gz=0'
		params += ' F="%s"' % fqfile1
		params += ' F2="%s"' % fqfile2
		params += ' filename="%s"' % infile
		if infile.endswith('.sam'):
			params += ' inputformat=sam'
		params += ' T="%s"' % path.join(tmpdir, infile + '.tmp')
		cmd = '{{args.biobambam_bamtofastq}} %s' % params
		runcmd (cmd)
	elif tool == 'bedtools':
		cmd = '{{args.bedtools}} bamtofastq %s -i "%s" -fq "%s" -fq2 "%s"' % (params, infile, fqfile1, fqfile2)
		runcmd (cmd)
	elif tool == 'samtools':
		cmd = '{{args.samtools}} fastq %s -t -1 "%s" -2 "%s" "%s"' % (params, fqfile1, fqfile2, infile)
		runcmd (cmd)
	elif tool == 'picard':
		mem = memtoJava({{ args.mem | quote }})
		cmd = '{{args.picard}} SamToFastq %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" F="%s" F2="%s"' % (mem, tmpdir, tmpdir, infile, fqfile1, fqfile2)
		runcmd (cmd)
		
	if {{args.gz | bool}}:
		runcmd ('gzip "%s"' % (fqfile1))
		runcmd ('gzip "%s"' % (fqfile2))
except:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)
"""

"""
@name:
	pBam2FastqSE
@description:
	Convert sam/bam files to single-end fastq files.
@input:
	`infile:file`: The sam/bam file. 
		- Sam files only available for biobambam, picard
@output:
	`fqfile:file`: The fastq file
@args:
	`tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
	`biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
	`bedtools`            : The path of bedtools. Default: bedtools
	`samtools`            : The path of samtools. Default: samtools
	`picard`              : The path of picard. Default: picard
	`mem`                 : The memory to be used by picard. Default: 8G
	`gz`                  : Whether gzip the output files. Default: True
	`params`:             : Other params for `tool`. Default: ''
	`tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`
@requires:
	[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
	[biobambam](https://github.com/gt1/biobambam2)
	[samtools](https://github.com/samtools/samtools)
	[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
"""
pBam2FastqSE                           = proc (desc = 'Convert bam files to single-end fastq files.')
pBam2FastqSE.input                     = "infile:file"
pBam2FastqSE.output                    = "fqfile:file:{{ infile | fn | fn | fn }}.fq{{args.gz | lambda x: '.gz' if x else ''}}"
pBam2FastqSE.args.tool                 = 'biobambam'
pBam2FastqSE.args.biobambam_bamtofastq = 'bamtofastq'
pBam2FastqSE.args.bedtools             = 'bedtools'
pBam2FastqSE.args.samtools             = 'samtools'
pBam2FastqSE.args.picard               = 'picard'
pBam2FastqSE.args.mem                  = '8G' # only for picard
pBam2FastqSE.args.gz                   = True
pBam2FastqSE.args.params               = ''
pBam2FastqSE.args.tmpdir               = __import__('tempfile').gettempdir() 
pBam2FastqSE.args._runcmd              = runcmd.python
pBam2FastqSE.args._memtoJava           = mem.toJava.python
pBam2FastqSE.lang                      = 'python'
pBam2FastqSE.script                    = """
from os import makedirs, path

tmpdir = path.join({{args.tmpdir | quote}}, "tmp.{{proc.id}}.{{proc.tag}}.{{proc.suffix}}.{{job.index}}")
if not path.exists(tmpdir):
	makedirs(tmpdir)

infile  = {{infile | quote}}
fqfile  = {{fqfile | quote}}
if {{args.gz | bool}}:
	fqfile = fqfile[:-3]

{{args._runcmd}}
{{args._memtoJava}}

params  = {{args.params | quote}}
tool    = {{args.tool | quote}}
try:
	if tool == 'biobambam':
		params += ' gz=0'
		params += ' S="%s"' % fqfile
		params += ' filename="%s"' % infile
		if infile.endswith('.sam'):
			params += ' inputformat=sam'
		params += ' T="%s"' % path.join(tmpdir, infile + '.tmp')
		cmd = '{{args.biobambam_bamtofastq}} %s' % params
		runcmd (cmd)
	elif tool == 'bedtools':
		cmd = '{{args.bedtools}} bamtofastq %s -i "%s" -fq "%s"' % (params, infile, fqfile)
		runcmd (cmd)
	elif tool == 'samtools':
		cmd = '{{args.samtools}} fastq %s -t -s "%s" "%s"' % (params, fqfile, infile)
		runcmd (cmd)
	elif tool == 'picard':
		mem = memtoJava({{ args.mem | quote }})
		cmd = '{{args.picard}} SamToFastq %s -Djava.io.tmpdir="%s" TMP_DIR="%s" I="%s" F="%s"' % (mem, tmpdir, tmpdir, infile, fqfile)
		runcmd (cmd)
		
	if {{args.gz | bool}}:
		runcmd ('gzip "%s"' % (fqfile))
except:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	from shutil import rmtree
	rmtree (tmpdir)
"""

"""
@name:
	pBam2Counts
"""
pBam2Counts = proc (desc = 'Extract read counts from RNA-seq bam files.')
pBam2Counts.input = 'infile:file'
pBam2Counts.output = 'outfile:file:{{infile | fn}}.counts'
pBam2Counts.args.tool = 'htseq'
pBam2Counts.args.htseq_count = 'htseq-count'
pBam2Counts.args.params = ''
pBam2Counts.args.refgene = ''
pBam2Counts.args._runcmd = runcmd.python
pBam2Counts.beforeCmd = checkArgsRefgene.bash
pBam2Counts.lang = 'python'
pBam2Counts.script = """
{{args._runcmd}}

tool = '{{args.tool}}'
if tool == 'htseq':
	samtype = ''
	if {{infile | quote}}.endswith('.bam'):
		samtype = '-f bam'
	cmd = '{{args.htseq-count}} {{args.params}} -r pos %s "{{infile}}" "{{args.refgene}}" > "{{outfile}}"' % (samtype)
	runcmd (cmd)
"""