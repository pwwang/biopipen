from shutil import move, rmtree
from os import path, symlink, remove
from sys import stderr
from pyppl import Box
from bioprocs.utils import shell, mem2, cmdargs
from bioprocs.utils.shell import runcmd2, RuncmdException, Shell

infile     = {{ i.infile | quote }}
inprefix   = {{ i.infile | fn | quote}}
outfile    = {{ o.outfile | quote }}
joboutdir  = {{job.outdir | quote}}
infmt      = {{ i.infile | ext | [1:] | quote}}
tmpdir     = {{ args.tmpdir | quote}}
tmpdir     = path.join (tmpdir, {{proc.id, fn(i.infile), str(job.index+1) | : '.'.join(a) | quote }})
steps      = {{ args.steps | repr }}
params     = {{ args.params | repr }}
argsmem    = {{ args.mem | repr}}
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
shell.TOOLS.update(dict(
	biobambam = biobambam,
	sambamba  = sambamba,
	samtools  = samtools,
	picard    = picard,
	elprep    = elprep
))

if steps.rmdup:
	steps.markdup = True
if not path.exists(tmpdir):
	shell.mkdir(tmpdir, p = True)
if steps.recal and tool != 'elprep':
	raise ValueError('Step "recal" is only enabled by "elprep", use pBamRecal for other tools.')

subshell = Shell(subcmd = True)

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

	Shell(dash = '', equal = '=').biobambam(**params).run()

def run_sambamba():
	global infile
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		subshell.sambamba.view(S = True, f = 'bam', o = outfile, t = nthread, _ = infile).run()
	else:
		bamfile = outfile
		if infmt == 'sam':
			bamfile = path.join(joboutdir, inprefix + '.s2b.bam')
			subshell.sambamba.view(S = True, f = 'bam', o = bamfile, t = nthread, _ = infile).run()
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
			subshell.sambamba.sort(**params).run()
			if infile != {{i.infile | quote}}:
				shell.rm(f = True, _ = infile)
			infile = bamfile
		if steps.markdup:
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			subshell.sambamba.markdup(r = steps.rmdup, t = nthread, tmpdir = tmpdir, _ = [infile, bamfile]).run()
			if infile != {{i.infile | quote}}:
				shell.rm(f = True, _ = infile)
			infile = bamfile
		if steps.index:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			else:
				subshell.sambamba.index(t = nthread, _ = [infile, infile + '.bai'])
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_samtools():
	global infile
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		subshell.samtools.view(b = True, o = outfile, O = 'bam', _ = infile).run()
	else:
		bamfile = outfile
		if steps.sort:
			mem = mem2(argsmem, 'M')
			bamfile = path.join(joboutdir, inprefix + '.sorted.bam')
			subshell.samtools.sort(
				m = mem + 'M', 
				n = sortby == 'queryname',
				o = bamfile,
				T = tmpdir,
				O = 'bam',
				_ = infile,
				**{'@': nthread}
			).run()
			if infile != {{i.infile | quote}}:
				shell.rm(infile, f = True)
			infile = bamfile
		if steps.markdup or steps.rmdup:
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			subshell.rmdup(infile, bamfile).run()
			if infile != {{i.infile | quote}}:
				shell.rm(infile, f = True)
			infile = bamfile
		if steps.index:
			subshell.samtools.index(bamfile, outfile + '.bai')
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_picard():
	global infile
	mem = mem2(argsmem, '-jdict')
	mem['-Djava.io.tmpdir'] = tmpdir
	shellpicard = Shell(subcmd = True, dash = '', equal = '=').picard(**mem)
	if not (steps.sort or steps.index or steps.markdup or steps.rmdup):
		shellpicard.SamFormatConverter(TMP_DIR = tmpdir, I = infile, O = outfile).run()
	else:
		bamfile = outfile
		if steps.sort:
			bamfile = path.join(joboutdir, inprefix + '.sorted.bam')
			shellpicard.ShortSam(TMP_DIR = tmpdir, I = infile, O = bamfile, SO = sortby).run()
			if infile != {{i.infile | quote}}:
				shell.rm(f = True, _ = infile)
			infile = bamfile
		if steps.markdup:
			mfile = "/dev/null"
			bamfile = path.join(joboutdir, inprefix + '.dedup.bam')
			shellpicard.MarkDuplicates(REMOVE_DUPLICATES = 'true' if steps.rmdup else 'false', TMP_DIR = tmpdir, I = infile, O = bamfile, M = mfile).run()
			if infile != {{i.infile | quote}}:
				shell.rm(f = True, _ = infile)
			infile = bamfile
		if steps.index:
			shellpicard.BuildBamIndex(TMP_DIR = tmpdir, I = infile, O = outfile + '.bai').run()
		if infile != outfile:
			if path.exists(infile + '.bai'):
				shell.mv(infile + '.bai', outfile + '.bai')
			shell.mv(infile, outfile)

def run_elprep():
	elpsh = Shell(subcmd = True, equal = ' ').elprep

	params['log-path']          = joboutdir
	params['nr-of-threads']     = nthread
	params['sorting-order']     = sortby if steps.sort else 'keep'
	params['mark-duplicates']   = steps.markdup
	params['remove-duplicates'] = steps.rmdup
	params['']                  = [infile, outfile]
	if steps.markdup:
		params['mark-optical-duplicates'] = path.join(joboutdir, inprefix + '.opticaldups.txt')
	if steps.recal:
		params['bqsr']           = path.join(joboutdir, inprefix + '.bqsr.txt')
		params['bqsr-reference'] = ref + '.elprep'
		if knownSites:
			params['known-sites'] = knownSites
	
	elpsh.filter(**params).run()
	if steps.index:
		subshell.samtools.index(outfile, outfile + '.bai')

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
