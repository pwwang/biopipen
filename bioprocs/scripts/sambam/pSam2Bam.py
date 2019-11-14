from shutil import move, rmtree
from os import path, symlink, remove
from sys import stderr
from pyppl import Box
from bioprocs.utils import mem2, shell2 as shell

infile     = {{ i.infile | quote }}
inprefix   = {{ i.infile | fn | quote}}
outfile    = {{ o.outfile | quote }}
joboutdir  = {{job.outdir | quote}}
infmt      = {{ i.infile | ext | [1:] | quote}}
tmpdir     = {{ args.tmpdir | quote}}
tmpdir     = path.join(tmpdir, {{proc.id, fn(i.infile), str(job.index+1) | @join: "_" | quote }})
steps      = {{ args.steps | repr }}
params     = {{ args.params | repr }}
argsmem    = {{ args.mem | repr }}
sortby     = {{ args.sortby | quote}}
nthread    = {{ args.nthread | repr}}
ref        = {{ args.ref | repr}}
knownSites = {{ args.knownSites | repr}}
biobambam  = {{ args.biobambam | quote }}
elprep     = {{ args.elprep | quote}}
sambamba   = {{ args.sambamba | quote}}
samtools   = {{ args.samtools | quote}}
picard     = {{ args.picard | quote}}
tool       = {{ args.tool | quote }}


shell.load_config(
	biobambam = biobambam,
	sambamba  = sambamba,
	samtools  = samtools,
	picard    = picard,
	elprep    = elprep
)

if steps.rmdup:
	steps.markdup = True
if not path.exists(tmpdir):
	shell.mkdir(tmpdir, p = True)
if steps.recal and tool != 'elprep':
	raise ValueError('Step "recal" is only enabled by "elprep", use pBamRecal for other tools.')


def run_biobambam():
	mem = mem2(argsmem, 'M')
	if steps.index:
		params.index         = 1
		params.indexfilename = outfile + '.bai'

	params.I              = infile
	params.O              = outfile
	params.SO             = sortby
	params.blockme        = mem
	params.tmpfile        = path.join(tmpdir, 'biobambam.tmp')
	params.inputformat    = infmt
	params.outfmt         = 'bam'
	params.inputthreads   = nthread
	params.outputthreads  = nthread
	params.markduplicates = int(steps.markdup)
	params.rmdup          = int(steps.rmdup)

	shell.fg.biobambam(**params)

def run_sambamba():
	global infile
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		shell.fg.sambamba.view(S = True, f = 'bam', o = outfile, t = nthread, _ = infile)
	else:
		bamfile = outfile
		if infmt == 'sam':
			bamfile = path.join(joboutdir, inprefix + '.s2b.bam')
			shell.fg.sambamba.view(S = True, f = 'bam', o = bamfile, t = nthread, _ = infile)
			infile = bamfile
		if steps.sort:
			if sortby == 'queryname':
				params.n = True
				params.N = True
			bamfile       = path.join(joboutdir, inprefix + '.sorted.bam')
			params.m      = argsmem
			params.tmpdir = tmpdir
			params.o      = bamfile
			params.t      = nthread
			params._      = infile
			shell.fg.sambamba.sort(**params)
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.markdup:
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			shell.fg.sambamba.markdup(r = steps.rmdup, t = nthread, tmpdir = tmpdir,
				_ = [infile, bamfile])
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.index:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			else:
				shell.fg.sambamba.index(t = nthread, _ = [infile, infile + '.bai'])
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_samtools():
	global infile
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		shell.fg.samtools.view(b = True, o = outfile, O = 'bam', _ = infile)
	else:
		bamfile = outfile
		if steps.sort:
			mem = mem2(argsmem, 'M')
			bamfile = path.join(joboutdir, inprefix + '.sorted.bam')
			shell.fg.samtools.sort(
				m     = mem + 'M',
				n     = sortby == 'queryname',
				o     = bamfile,
				T     = tmpdir,
				O     = 'bam',
				_     = infile,
				**{'@': nthread}
			)
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.markdup or steps.rmdup:
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			shell.fg.samtools.rmdup(infile, bamfile)
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.index:
			shell.samtools.index(bamfile, outfile + '.bai')
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_picard():
	global infile
	mem = mem2(argsmem, '-jdict')
	mem['-Djava.io.tmpdir'] = tmpdir
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		shell.fg.picard.SamFormatConverter(TMP_DIR = tmpdir, I = infile, O = outfile)
	else:
		bamfile = outfile
		if steps.sort:
			bamfile = path.join(joboutdir, inprefix + '.sorted.bam')
			shell.fg.picard.ShortSam(TMP_DIR = tmpdir, I = infile, O = bamfile, SO = sortby)
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.markdup:
			mfile = "/dev/null"
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			shell.fg.picard.MarkDuplicates(REMOVE_DUPLICATES = 'true' if steps.rmdup else 'false',
				TMP_DIR = tmpdir, I = infile, O = bamfile, M = mfile)
			if infile != {{i.infile | quote}}:
				shell.rm_rf(infile)
			infile = bamfile
		if steps.index:
			shell.fg.picard.BuildBamIndex(TMP_DIR = tmpdir, I = infile, O = outfile + '.bai')
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_elprep():

	params['log-path']                       = joboutdir
	params['nr-of-threads']                  = nthread
	params['sorting-order']                  = sortby if steps.sort else 'keep'
	params['mark-duplicates']                = steps.markdup
	params['remove-duplicates']              = steps.rmdup
	params['intermediate-files-output-type'] = 'bam'
	params['contig-group-size']              = 1
	params['tmp-path']                       = tmpdir
	params['']                               = [infile, outfile]
	if steps.markdup:
		params['mark-optical-duplicates'] = path.join(joboutdir, inprefix + '.opticaldups.txt')
	if steps.recal:
		params['bqsr']           = path.join(joboutdir, inprefix + '.bqsr.txt')
		params['bqsr-reference'] = ref + '.elprep'
		if knownSites:
			params['known-sites'] = knownSites

	shell.fg.elprep.sfm(**params)
	if steps.index:
		shell.samtools.index(outfile, outfile + '.bai')

tools = {
	'biobambam': run_biobambam,
	'sambamba' : run_sambamba,
	'samtools' : run_samtools,
	'picard'   : run_picard,
	'elprep'   : run_elprep
}

try:
	tools[tool]()
except KeyError:
	raise KeyError('Tool {!r} not supported.'.format(tool))
except Exception as ex:
	stderr.write ("Job failed: %s" % str(ex))
	raise
finally:
	rmtree (tmpdir)
