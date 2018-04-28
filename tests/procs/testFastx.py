import unittest, testly, helpers
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config, testdirs
from bioprocs.fastx import pFastqSim, pFastQC, pFastMC, pFastqTrim, pFastqSETrim, pFastqSE2Sam, pFastq2Sam, pFastq2Expr
from bioprocs.common import pFiles2Dir
from bioprocs import params

class TestFastx (helpers.TestCase):

	testdir, indir, outdir = testdirs('TestFastx')

	def dataProvider_test1_FastqPESim(self):
		yield 't1', 'wgsim', path.join(self.outdir, 'fastqsim_wgsim_1.fastq')
		yield 't2', 'dwgsim',  path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')

	def test1_FastqPESim (self, tag, tool, outfile):
		pFastqSimTest = pFastqSim.copy(tag = tag)
		pFastqSimTest.input = [8525]
		pFastqSimTest.args.ref = params.ref.value
		pFastqSimTest.args.num   = 100
		pFastqSimTest.args.gz    = False
		pFastqSimTest.args.tool  = tool
		PyPPL(config).start(pFastqSimTest).run()
		self.assertFileEqual(pFastqSimTest.channel.get(), outfile)

	def dataProvider_test2_FastQC(self):
		infile = path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')
		yield 't1', infile
		infile = path.join(self.outdir, 'fastqsim_wgsim_1.fastq')
		yield 't2', infile

	def test2_FastQC (self, tag, infile):
		pFastQCTest = pFastQC.copy(tag = tag)
		pFastQCTest.input = [infile]
		PyPPL(config).start(pFastQCTest).run()
		self.assertTrue(path.exists(path.join(pFastQCTest.channel.get(), path.splitext(path.basename(infile))[0] + '_fastqc.html')))

	def dataProvider_test2_pFastq2Sam(self):
		infile1 = path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')
		infile2 = path.join(self.outdir, 'fastqsim_dwgsim_2.fastq')
		yield 't1', infile1, infile2, 'bwa', path.join(self.outdir, 'fastq2sam_bwa.bam')
		yield 't2', infile1, infile2, 'ngm', path.join(self.outdir, 'fastq2sam_ngm.bam')
		yield 't3', infile1, infile2, 'bowtie2', path.join(self.outdir, 'fastq2sam_bowtie2.bam')
		yield 't4', infile1, infile2, 'star', path.join(self.outdir, 'fastq2sam_star.bam')

	def test2_pFastq2Sam (self, tag, infile1, infile2, tool, outfile):
		pFastq2SamTest = pFastq2Sam.copy(tag = tag)
		pFastq2SamTest.args.nthread = 2
		pFastq2SamTest.input        = (infile1, infile2)
		pFastq2SamTest.args.outfmt  = 'bam'
		pFastq2SamTest.args.ref     = params.ref.value
		pFastq2SamTest.args.tool    = tool
		PyPPL(config).start(pFastq2SamTest).run()
		self.assertBamCountEqual(pFastq2SamTest.channel.get(), outfile)

	def dataProvider_test2_pFastqSE2Sam(self):
		infile = path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')
		yield 't1', infile, 'bwa', path.join(self.outdir, 'fastqse2sam_bwa.bam')
		yield 't2', infile, 'ngm', path.join(self.outdir, 'fastqse2sam_ngm.bam')
		yield 't3', infile, 'bowtie2', path.join(self.outdir, 'fastqse2sam_bowtie2.bam')
		yield 't4', infile, 'star', path.join(self.outdir, 'fastqse2sam_star.bam')

	def test2_pFastqSE2Sam (self, tag, infile, tool, outfile):
		pFastqSE2SamTest = pFastqSE2Sam.copy(tag = tag)
		pFastqSE2SamTest.args.nthread = 2
		pFastqSE2SamTest.input        = [infile]
		pFastqSE2SamTest.args.outfmt  = 'bam'
		pFastqSE2SamTest.args.ref     = params.ref.value
		pFastqSE2SamTest.args.tool    = tool
		PyPPL(config).start(pFastqSE2SamTest).run()
		self.assertBamCountEqual(pFastqSE2SamTest.channel.get(), outfile)

	def dataProvider_test3_MultiQC(self):
		infiles = [
			path.join(self.outdir, 'fastqc-out1'),
			path.join(self.outdir, 'fastqc-out2')
		]
		outfile = path.join(self.outdir, 'multiqc-out')
		yield 't1', infiles, outfile

	def test3_MultiQC (self, tag, infiles, outfile):
		pFiles2DirFastMC = pFiles2Dir.copy()
		pFiles2DirFastMC.input = [infiles]
		pFastMCTest = pFastMC.copy(tag = tag)
		pFastMCTest.depends = pFiles2DirFastMC
		PyPPL(config).start(pFiles2DirFastMC).run()

	def dataProvider_test4_PETrim(self):
		infiles = [
			path.join(self.outdir, 'fastqsim_dwgsim_1.fastq'),
			path.join(self.outdir, 'fastqsim_dwgsim_2.fastq'),
		]
		yield infiles[0], infiles[1], 'trimmomatic', path.join(self.outdir, 'fastqtrim-pe-trimmometic_1.fastq.gz')
		yield infiles[0], infiles[1], 'cutadapt', path.join(self.outdir, 'fastqtrim-pe-cutadapt_1.fastq.gz')
		yield infiles[0], infiles[1], 'skewer', path.join(self.outdir, 'fastqtrim-pe-skewer_1.fastq.gz')

	def test4_PETrim (self, infile1, infile2, tool, outfile):
		pFastqTrimTest = pFastqTrim.copy(tag = tool)
		pFastqTrimTest.input = (infile1, infile2)
		pFastqTrimTest.args.gz = True
		pFastqTrimTest.args.tool = tool
		pFastqTrimTest.args.nthread = 2

		PyPPL(config).start(pFastqTrimTest).run()
		self.assertFileEqual(pFastqTrimTest.channel.get(), outfile)

	def dataProvider_test4_SETrim(self):
		infile = path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')
		yield infile, 'trimmomatic', path.join(self.outdir, 'fastqtrim-se-trimmometic_1.fastq.gz')
		yield infile, 'cutadapt', path.join(self.outdir, 'fastqtrim-se-cutadapt_1.fastq.gz')
		yield infile, 'skewer', path.join(self.outdir, 'fastqtrim-se-skewer_1.fastq.gz')

	def test4_SETrim (self, infile, tool, outfile):
		pFastqSETrimTest = pFastqSETrim.copy(tag = tool)
		pFastqSETrimTest.input = [infile]
		pFastqSETrimTest.args.gz = True
		pFastqSETrimTest.args.nthread = 2
		pFastqSETrimTest.args.tool = tool

		PyPPL(config).start(pFastqSETrimTest).run()
		self.assertFileEqual(pFastqSETrimTest.channel.get(), outfile)

	def dataProvider_test5_pFastq2Expr(self):
		infile1 = path.join(self.outdir, 'fastqsim_dwgsim_1.fastq')
		infile2 = path.join(self.outdir, 'fastqsim_dwgsim_2.fastq')
		outfile = path.join(self.outdir, 'fastq2expr.txt')
		yield 't1', infile1, infile2, outfile

	def test5_pFastq2Expr(self, tag, infile1, infile2, outfile):
		pFastq2ExprTest = pFastq2Expr.copy(tag = tag)
		pFastq2ExprTest.input = infile1, infile2
		pFastq2ExprTest.args.nthread = 10
		PyPPL(config).start(pFastq2ExprTest).run()
		self.assertFileEqual(pFastq2ExprTest.channel.outfile.get(), outfile)

if __name__ == '__main__':
	testly.main(failfast=True)
