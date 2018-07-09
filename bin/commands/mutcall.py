#!/usr/bin/env python
# Call mutations from sequencing data.
import sys
from os import path
from pyppl import PyPPL, Channel, Box
from bioprocs import params
from bioprocs.common import pFiles2Dir
from bioprocs.sambam import pBam2Gmut, pBamPair2Smut
from bioprocs.utils.tsvio import TsvReader
from bioprocs.utils.sampleinfo import SampleInfo
from bioaggrs.wxs import aEBam2Bam, aFastq2Bam, aFastqSE2Bam, aBam2SCnv, aBam2GCnv

#params.prefix('-')
params.intype         = 'bam' # ebam, fastq
params.intype.desc    = 'The input file types. Either bam, ebam or fastq.\nEbam means the bam files need to reformatted into fastq files.'
params.muts           = ['germ'] # or ['germ', 'soma', 'scnv', 'gcnv']
params.muts.type      = list
params.muts.desc      = 'What kind of mutations to call.\nNote: soma need paired information'
params.indir.required = True
params.indir.desc     = 'The input directory containing input files.'
params.saminfo.required = True
params.saminfo.desc = """The sample information file:
Column 1: the basename of the sample file in '-indir'
Column 2: Group if '-muts' includes 'soma', basically, 'TUMOR' and 'NORMAL'
Column 3: Patient if multiple samples belong to the same patient.
Example:
Sample	Group	Patient
A1.bam	TUMOR	A
A2.bam	NORMAL	A
B1.bam	TUMOR	A
B2.bam	NORMAL	A
"""

params.exdir.required = True
params.exdir.desc     = 'Where to export the result files.'
params.runner         = 'local'
params.runner.desc    = 'The runner to run the processes.'
params.forks          = 1
params.forks.desc     = 'How many jobs to run simultanuously.'
params.logfile        = ''
params.logfile.desc   = 'Where to save the logs.'
params.compress       = True
params.compress.desc  = 'Use gzip and bam file to save space.'

def _procconfig(kwargs = None):
	params = _getparams(kwargs or {})

	starts = []
	
	saminfo = SampleInfo(params.saminfo)

	aEBam2Bam.pFastq2Sam.args.tool    = 'bowtie2'

	pBamDir         = pFiles2Dir
	pBamDir.runner  = 'local'
	if params.intype == 'ebam':
		#aEBam2Bam.input = [Channel.fromPattern(path.join(params.indir, '*.bam'))]
		aEBam2Bam.input = [Channel.create(saminfo.toChannel(params.indir)).unique()]
		if params.compress:
			aEBam2Bam.args.gz = True
			aEBam2Bam.pFastq2Sam.args.outfmt = 'bam'

		pBamDir.depends = aEBam2Bam
		pBamDir.input   = lambda ch: [ch.flatten()]

		starts.append(aEBam2Bam)

	elif params.intype == 'sfq' or params.intype == 'sfastq':
		# single-end fastq files
		# *.fq, *.fq.gz *.fastq, *.fastq.gz
		# sample info should be:
		# +--------------+----------+---------+
		# | Sample	     | Patient  | Group   |
		# | x_Tumor.bam	 | x        | TUMOR	  |
		# | x_Normal.bam | x        | NORMAL  |
		# | ...          | ...      | ...     |
		# +--------------+----------+---------+
		# corresponding fastq files would be:
		# x_Tumor.fq  / x_Tumor.fastq  / x_Tumor.fq.gz  / x_Tumor.fastq.gz
		# x_Normal.fq / x_Normal.fastq / x_Normal.fq.gz / x_Normal.fastq.gz
		def bam2fq(fq):
			fqdir   = path.dirname(fq)
			bname   = path.splitext(path.basename(fq))[0]
			exts    = ['.fq', '.fq.gz', '.fq.gz', '.fastq.gz']
			fqfiles = [path.join(fqdir, bname + ext) for ext in exts]
			return [fqfile for fqfile in fqfiles if path.isfile(fqfile)][0]
			
		aFastqSE2Bam.input = [bam2fq(fq) for fq in Channel.create(saminfo.toChannel(params.indir)).unique()]

		pBamDir.depends = aFastqSE2Bam.pBamRecal
		pBamDir.input   = lambda ch: [ch.flatten()]

		starts.append(aFastqSE2Bam)
	elif params.intype == 'fq' or params.intype == 'fastq':
		# pair-end fastq files
		# *.fq, *.fq.gz *.fastq, *.fastq.gz
		# sample info should be:
		# +--------------+----------+---------+
		# | Sample	     | Patient  | Group   |
		# | x_Tumor.bam	 | x        | TUMOR	  |
		# | x_Normal.bam | x        | NORMAL  |
		# | ...          | ...      | ...     |
		# +--------------+----------+---------+
		# corresponding fastq files would be:
		# x_Tumor_1.fq  / x_Tumor_1.fastq  / x_Tumor_1.fq.gz  / x_Tumor_1.fastq.gz
		# x_Tumor_2.fq  / x_Tumor_2.fastq  / x_Tumor_2.fq.gz  / x_Tumor_2.fastq.gz
		# x_Normal_1.fq / x_Normal_1.fastq / x_Normal_1.fq.gz / x_Normal_1.fastq.gz
		# x_Normal_2.fq / x_Normal_2.fastq / x_Normal_2.fq.gz / x_Normal_2.fastq.gz
		def bam2fqpair(fq):
			fqdir    = path.dirname(fq)
			bname    = path.splitext(path.basename(fq))[0]
			exts1    = ['_1.fq', '_1.fq.gz', '_1.fq.gz', '_1.fastq.gz']
			exts2    = ['_2.fq', '_2.fq.gz', '_2.fq.gz', '_2.fastq.gz']
			fqfiles1 = [path.join(fqdir, bname + ext) for ext in exts1]
			fqfile1  = [fqfile for fqfile in fqfiles1 if path.isfile(fqfile)][0]
			fqfiles2 = [path.join(fqdir, bname + ext) for ext in exts2]
			fqfile2  = [fqfile for fqfile in fqfiles2 if path.isfile(fqfile)][0]
			return fqfile1, fqfile2
			
		aFastqSE2Bam.input = [bam2fqpair(fq) for fq in saminfo.toChannel(params.indir)]

		pBamDir.depends = aFastqSE2Bam.pBamRecal
		pBamDir.input   = lambda ch: [ch.flatten()]

		starts.append(aFastqSE2Bam)
	else:
		pBamDir.input = [saminfo.toChannel(params.indir)]
		starts.append(pBamDir)

	if 'germ' in params.muts:
		pBam2Gmut.depends = pBamDir
		pBam2Gmut.input   = lambda ch: ch.expand(0, "*.bam")
		pBam2Gmut.exdir   = path.join(params.exdir, 'germline')
	if 'soma' in params.muts:
		pBamPair2Smut.depends = pBamDir
		pBamPair2Smut.input   = lambda ch: saminfo.toChannel(ch.get(), paired = True)
		pBamPair2Smut.exdir   = path.join(params.exdir, 'somatic')
	if 'scnv' in params.muts:
		aBam2SCnv.pBamDir.depends   = pBamDir
		aBam2SCnv.pSampleInfo.input = [params.saminfo]
		aBam2SCnv.exdir             = path.join(params.exdir, 'scnv')
		starts.append(aBam2SCnv)
	if 'gcnv' in params.muts:
		aBam2GCnv.pBamDir.depends   = pBamDir
		aBam2GCnv.pSampleInfo.input = [params.saminfo]
		aBam2GCnv.exdir             = path.join(params.exdir, 'gcnv')
		starts.append(aBam2GCnv)
		
	config = {
		'default': {'forks': int(params.forks)},
		'_log' : {'file': params.logfile}
	}
	return starts, config, params.runner

def _getparams(kwargs):
	global params
	
	if len(sys.argv) > 1 and sys.argv[1] == path.splitext(path.basename(__file__))[0]:
		# called from api
		params('hopts', '', True)
		for key, val in kwargs.items():
			setattr(params, key, val)
		return params.parse(args = []).asDict()
	else:
		# called directly
		return params.parse().asDict()

def run(*args, **kwargs):
	proc, config, runner = _procconfig(kwargs)
	PyPPL(config).start(proc).run(runner)

if __name__ == '__main__':
	run()

