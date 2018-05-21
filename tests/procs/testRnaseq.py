import unittest, helpers, testly
from os import path
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config
from bioprocs.rnaseq import pEXPRdir2Matrix, pBatchEffect, pRawCounts2, pRNAseqDEG, pCoexp, p2RawCounts

class TestRnaseq (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestRnaseq')

	def test1pEXPRdir2Matrix(self):
		pEXPRdir2Matrix.input = [getfile()]
		pEXPRdir2Matrix.expect = "grep Gene100 {{out.outfile}} && grep Sample20 {{out.outfile}} && !(grep Composite {{out.outfile}})"
		pEXPRdir2Matrix.args.pattern  = '*.expr'
		pEXPRdir2Matrix.args.plot.boxplot  = True
		pEXPRdir2Matrix.args.plot.heatmap  = True
		pEXPRdir2Matrix.args.plot.histogram = True
		PyPPL().start(pEXPRdir2Matrix).run()
		procOK(pEXPRdir2Matrix, 'expdir2mat.txt', self)

	def test2pBatchEffect(self):
		pBatchEffect.input = [(getfile('expdir2mat.txt', input=False), getfile('batch.txt'))]
		pBatchEffect.args.plot.boxplot  = True
		pBatchEffect.args.plot.heatmap  = True
		pBatchEffect.args.plot.histogram = True
		PyPPL().start(pBatchEffect).run()
		procOK(pBatchEffect, 'batcheffect.txt', self)

	def test3pRawCounts2(self):
		pRawCounts2.input = [getfile('expdir2mat.txt', input=False)]
		pRawCounts2.args.plot.boxplot  = True
		pRawCounts2.args.plot.heatmap  = True
		pRawCounts2.args.plot.histogram = True
		PyPPL().start(pRawCounts2).run()
		procOK(pRawCounts2, 'rawcounts2.txt', self)

	def dataProvider_test4pRNAseqDEG(self):
		infile   = path.join(self.outdir, 'expdir2mat.txt')
		samfile  = path.join(self.indir, 'batch.txt')

		pPos = Proc(desc = 'Make values positive.')
		pPos.input  = {"infile:file": infile}
		pPos.output = "outfile:file:{{in.infile|bn}}"
		pPos.script = "sed 's/-//g' {{in.infile}} > {{out.outfile}}"

		pPairfile = Proc(desc = 'Add patient to group file.')
		pPairfile.input  = {"infile:file": samfile}
		pPairfile.output = "outfile:file:{{in.infile|bn}}"
		pPairfile.script = 'paste {{in.infile}} <(echo -e "Patient\\n$(seq 1 10)\\n$(seq 1 10)") > {{out.outfile}}'

		plot     = {'maplot': True, 'mdsplot':True, 'volplot': True, 'heatmap': True}

		tag1     = 'edger'
		tool1    = 'edger'
		outfile1 = path.join(self.outdir, 'deg-edger.txt')

		yield tag1, pPos, pPairfile, plot, tool1, outfile1, False, samfile

		tag2     = 'deseq2'
		tool2    = 'deseq2'
		outfile2 = path.join(self.outdir, 'deg-deseq2.txt')
		yield tag2, pPos, pPairfile, plot, tool2, outfile2, False, samfile

		tag3     = 'edgerp'
		tool3    = 'edger'
		outfile3 = path.join(self.outdir, 'deg-edger-paired.txt')
		yield tag3, pPos, pPairfile, plot, tool3, outfile3, True, samfile

		tag4     = 'deseq2p'
		tool4    = 'deseq2'
		outfile4 = path.join(self.outdir, 'deg-deseq2-paired.txt')
		yield tag4, pPos, pPairfile, plot, tool4, outfile4, True, samfile

	def test4pRNAseqDEG(self, tag, pPos, pPairfile, plot, tool, outfile, paired, samfile):
		pPos = pPos.copy(tag = tag)
		pPairfile = pPairfile.copy(tag = tag)
		
		pRNAseqDEGTest = pRNAseqDEG.copy(tag = tag)
		pRNAseqDEGTest.depends   = (pPos, pPairfile) if paired else pPos
		pRNAseqDEGTest.args.plot = plot
		pRNAseqDEGTest.args.tool = tool
		if not paired:
			pRNAseqDEGTest.input = lambda ch: ch.cbind(samfile)

		PyPPL(config).start(pPos, pPairfile).run()
		self.assertFileEqual(pRNAseqDEGTest.channel.outfile.get(), outfile)

	def testCoexp(self):
		pCoexp.input = [getfile('coexp.expr.txt')]
		pCoexp.args.pval = True
		PyPPL().start(pCoexp).run()
		procOK(pCoexp, 'coexp.txt', self)

	def test2RawCounts(self):
		p2RawCounts.input = [getfile('p2rawcounts.expr.txt')]
		p2RawCounts.args.refgene = getfile('p2rawcounts.glen.txt')
		PyPPL().start(p2RawCounts).run()
		procOK(p2RawCounts, 'p2rawcounts.txt', self)

if __name__ == '__main__':
	testly.main(failfast=True)
