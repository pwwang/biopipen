import unittest
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config
from bioprocs.rnaseq import pExpdir2Matrix, pBatchEffect, pRawCounts2, pRnaseqDeg, pCoexp, p2RawCounts

class TestRnaseq (unittest.TestCase):

	def test1pExpdir2Matrix(self):
		pExpdir2Matrix.input = [getfile()]
		pExpdir2Matrix.expect = "grep Gene100 {{out.outfile}} && grep Sample20 {{out.outfile}} && !(grep Composite {{out.outfile}})"
		pExpdir2Matrix.args.pattern  = '*.expr'
		pExpdir2Matrix.args.boxplot  = True
		pExpdir2Matrix.args.heatmap  = True
		pExpdir2Matrix.args.histplot = True
		PyPPL().start(pExpdir2Matrix).run()
		procOK(pExpdir2Matrix, 'expdir2mat.txt', self)

	def test2pBatchEffect(self):
		pBatchEffect.input = [(getfile('expdir2mat.txt', input=False), getfile('batch.txt'))]
		pBatchEffect.args.boxplot  = True
		pBatchEffect.args.heatmap  = True
		pBatchEffect.args.histplot = True
		PyPPL().start(pBatchEffect).run()
		procOK(pBatchEffect, 'batcheffect.txt', self)

	def test3pRawCounts2(self):
		pRawCounts2.input = [getfile('expdir2mat.txt', input=False)]
		pRawCounts2.args.boxplot  = True
		pRawCounts2.args.heatmap  = True
		pRawCounts2.args.histplot = True
		PyPPL().start(pRawCounts2).run()
		procOK(pRawCounts2, 'rawcounts2.txt', self)

	def test4pRnaseqDeg(self):
		pPos = Proc(desc = 'Make values positive.')
		pPos.input  = {"infile:file": getfile('expdir2mat.txt', input=False)}
		pPos.output = "outfile:file:{{in.infile|fn}}"
		pPos.script = "sed 's/-//g' {{in.infile}} > {{out.outfile}}"

		pPairfile = Proc(desc = 'Add patient to group file.')
		pPairfile.input = {"infile:file": getfile('batch.txt')}
		pPairfile.output = "outfile:file:{{in.infile|fn}}"
		pPairfile.script = 'paste {{in.infile}} <(echo -e "Patient\\n$(seq 1 10)\\n$(seq 1 10)") > {{out.outfile}}'

		pRnaseqDegDeseq         = pRnaseqDeg.copy()
		pRnaseqDegPaired        = pRnaseqDeg.copy()
		pRnaseqDegDSPaired      = pRnaseqDeg.copy()
		pRnaseqDeg.depends      = pPos
		pRnaseqDeg.args.heatmap = True
		pRnaseqDeg.args.maplot  = True
		pRnaseqDeg.input        = lambda ch: ch.cbind([getfile('batch.txt')])

		pRnaseqDegDeseq.depends      = pPos
		pRnaseqDegDeseq.args.heatmap = True
		pRnaseqDegDeseq.args.tool    = 'deseq2'
		pRnaseqDegDeseq.input        = lambda ch: ch.cbind([getfile('batch.txt')])

		pRnaseqDegPaired.depends     = pPos, pPairfile
		pRnaseqDegDSPaired.depends   = pPos, pPairfile
		pRnaseqDegDSPaired.args.tool = 'deseq2'
		PyPPL().start(pPos, pPairfile).run()
		procOK(pRnaseqDeg, 'deg-edger.txt', self)
		procOK(pRnaseqDegDeseq, 'deg-deseq2.txt', self)
		procOK(pRnaseqDegPaired, 'deg-edger-paired.txt', self)
		procOK(pRnaseqDegDSPaired, 'deg-deseq2-paired.txt', self)

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
	unittest.main(failfast=True)