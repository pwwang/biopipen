import re
from os import path, symlink
from sys import stdout, stderr
from time import sleep
from shutil import move
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs, reference, logger
from bioprocs.utils.poll import Poll

# variables
infile   = {{i.fq | quote}}
outfile  = {{o.outfile | quote}}
outfmt   = {{args.outfmt | quote}}
outdir   = {{job.outdir | quote}}
ref      = {{args.ref | quote}}
refgene  = {{args.refgene | quote}}
ref2     = path.join(outdir, path.basename(ref))
workdir  = {{proc.workdir | quote}}
rg       = {{args.rg}}
params   = {{args.params}}
samtools = {{args.samtools | quote}}
bowtie2  = {{args.bowtie2 | quote}}
bwa      = {{args.bwa | quote}}
ngm      = {{args.ngm | quote}}
star     = {{args.star | quote}}
nthread  = {{args.nthread}}
joblen   = {{proc.size}}
jobindex = {{job.index}}
bowtie2_build  = {{args.bowtie2_build | quote}}
poll     = Poll(workdir, joblen, jobindex)

# determine whether the reference file is specified
logger.info('Checking if reference exists: %s' % ref)
reference.check(ref)

# detemine default read group
rg = {key.upper():val for key, val in rg.items()}
if not rg['ID']:
	g = re.search (r'[^a-zA-Z0-9]+(L\\d+)[^a-zA-Z0-9]+', {{o.outfile | fn | quote}})
	rg['ID'] = g.group(1) if g else "{{o.outfile | fn}}.L{{job.index}}"
if not rg['SM']:
	rg['SM'] = {{o.outfile | fn | quote}}

def sam2bam(samfile, bamfile):
	logger.info('Converting sam to bam: ')
	logger.info('- %s' % samfile)
	logger.info('- %s' % bamfile)
	cmd = '%s view -Sb "%s" > "%s"; rm -f "%s"' % (samtools, samfile, bamfile, samfile)
	runcmd(cmd)

def checkAndBuildIndex(indexes, buildcmd, buildcmd2):
	def todo(indexes, ref, buildcmd, ref2, buildcmd2):
		logger.info('Checking if reference file is indexed.')
		if not reference.checkIndex(indexes):
			logger.info('No, trying to build index at the first job.')
			reference.buildIndex(ref, buildcmd, ref2, buildcmd2)
	poll.first(todo, indexes, ref, buildcmd, ref2, buildcmd2)
	return ref2 if path.exists(ref2) else ref

############ bowtie2
{% case args.tool %}
{% when 'bowtie2' %}
# check reference index files
indexes = [
	"%s.1.bt2" % ref,
	"%s.2.bt2" % ref,
	"%s.3.bt2" % ref,
	"%s.4.bt2" % ref,
	"%s.rev.1.bt2" % ref,
	"%s.rev.2.bt2" % ref
]
buildcmdbase = "{bowtie2_build} --thread {nthread} '{ref}' '{ref}'"
buildcmd     = buildcmdbase.format(
	bowtie2_build = bowtie2_build, nthread = nthread, ref = ref
)
buildcmd2    = buildcmdbase.format(
	bowtie2_build = bowtie2_build, nthread = nthread, ref = ref2
)
ref = checkAndBuildIndex(indexes, buildcmd, buildcmd2)

# do mapping
params['threads'] = nthread
params['rg-id']   = rg['ID']
i = 0
for k,v in rg.items():
	if k == 'ID': continue
	params['rg' + ' '*i] = k + ':' + v
	i += 1
params['x'] = ref
params['U'] = infile
params['S'] = path.splitext(outfile)[0] + '.sam'

cmd = '{bowtie2} {args}'.format(bowtie2 = bowtie2, args = cmdargs(params, equal = ' '))
logger.info('Do mapping ...')
runcmd (cmd)
if outfmt == 'bam':	sam2bam(params['S'], outfile)

############ bwa
{% when 'bwa' %}
indexes = [
	"%s.bwt" % ref,
	"%s.amb" % ref,
	"%s.ann" % ref,
	"%s.pac" % ref,
	"%s.sa" % ref,
]
buildcmdbase = "{bwa} index '{ref}'"
buildcmd     = buildcmdbase.format(bwa = bwa, ref = ref)
buildcmd2    = buildcmdbase.format(bwa = bwa, ref = ref2)
ref = checkAndBuildIndex(indexes, buildcmd, buildcmd2)

# do mapping
params['t'] = nthread
params['R'] = "@RG\\tID:%s\\t%s" % (rg['ID'], "\\t".join(k + ':' +v for k,v in rg.items() if k!='ID'))
outfile2    = path.splitext(outfile)[0] + '.sam'
cmd         = "{bwa} mem {args} '{ref}' '{infile}' > '{outfile}'".format(
	bwa     = bwa,
	args    = cmdargs(params),
	ref     = ref,
	infile  = infile,
	outfile = outfile2
)
logger.info('Do mapping ...')
runcmd (cmd)
if outfmt == 'bam':	sam2bam(outfile2, outfile)

{% when 'ngm' %}
indexes = [
	"%s-enc.2.ngm" % ref,
	"%s-ht-13-2.3.ngm" % ref
]
buildcmdbase = "{ngm} -r '{ref}' -t {nthread}"
buildcmd     = buildcmdbase.format(ngm = ngm, ref = ref, nthread = nthread)
buildcmd2    = buildcmdbase.format(ngm = ngm, ref = ref2, nthread = nthread)
ref = checkAndBuildIndex(indexes, buildcmd, buildcmd2)

# do mapping
params['q'] = infile
params['r'] = ref
params['o'] = outfile
params['b'] = {{ args.outfmt | lambda x: x=='bam' }}
params.update({('rg-' + k.lower()):v for k,v in rg.items()})
params['t'] = nthread
cmd = '%s %s' % (ngm, cmdargs(params, equal = ' '))
logger.info('Do mapping ...')
runcmd (cmd)

############ star
{% when 'star' %}
indexes = [path.join(path.splitext(ref)[0] + '.star', x) for x in [
	"chrLength.txt",
	"chrNameLength.txt",
	"chrName.txt",
	"chrStart.txt",
	"exonGeTrInfo.tab",
	"exonInfo.tab",
	"geneInfo.tab",
	"Genome",
	"genomeParameters.txt",
	"SA",
	"SAindex",
	"sjdbInfo.txt",
	"sjdbList.fromGTF.out.tab",
	"sjdbList.out.tab",
	"transcriptInfo.tab"
]]
buildcmdbase = "mkdir -p '{stardir}'; {star} --runMode genomeGenerate --genomeDir '{stardir}' --genomeFastaFiles " \
               + "'{ref}' --sjdbOverhang 100 --sjdbGTFfile '{refgene}' --runThreadN {nthread}"
buildcmd = buildcmdbase.format(
	stardir = path.splitext(ref)[0] + '.star',
	star    = star,
	ref     = ref,
	refgene = refgene,
	nthread = nthread
)
buildcmd2 = buildcmdbase.format(
	stardir = path.splitext(ref2)[0] + '.star',
	star    = star,
	ref     = ref2,
	refgene = refgene,
	nthread = nthread
)
ref = checkAndBuildIndex(indexes, buildcmd, buildcmd2)
refDir = path.splitext(ref)[0] + '.star'
rfcmd  = "zcat" if infile.endswith('.gz') else 'bzcat' if infile.endswith('.bz2') else 'cat'

params['genomeDir']         = refDir
params['readFilesIn']       = infile
params['readFilesCommand']  = rfcmd
params['readNameSeparator'] = '.'
params['outFileNamePrefix'] = outdir + '/'
params['outSAMtype']        = [outfmt.upper(), 'Unsorted']
cmd                         = '%s %s' % (star, cmdargs(params, equal = ' '))
logger.info('Do mapping ...')
runcmd (cmd)

outfile2 = "{outdir}/Aligned.out.{outfmt}".format(outdir = outdir, outfmt = outfmt)
if path.exists(outfile2): move(outfile2, outfile)
{% endcase %}
