import helpers, testly
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.stats import pMetaPval, pMetaPval1, pSurvival, pChiSquare, pFisherExact, pPWFisherExact, pMediation, pHypergeom
from bioprocs.common import pFiles2Dir
from bioprocs.tsv import pRbind

class TestStats (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestStats')
	
	def dataProvider_test1MetaPval(self):
		yield 't1', self.indir, path.join(self.outdir, 'metapval.txt'), {'pattern': "testMetaPval*.txt"}

	def test1MetaPval(self, tag, infile, outfile, args = None):
		pMetaPvalTest = pMetaPval.copy(tag = tag)
		pMetaPvalTest.input = [infile]
		pMetaPvalTest.args.update(args or {})
		PyPPL(config).start(pMetaPvalTest).run()
		self.assertFileEqual(pMetaPvalTest.channel.get(), outfile)

	def test2MetaPval1(self):
		pRbind.input = [[getfile('testMetaPval%s.txt' % i) for i in range(10)]]
		pRbind.args.inopts.rnames = False
		pRbind.args.inopts.cnames = True

		pMetaPval1.depends = pRbind
		
		PyPPL(config).start(pRbind).run()
		procOK(pMetaPval1, 'metapval2.txt', self)

	def dataProvider_testSurvival(self):
		infile1  = path.join(self.indir, 'survival1.txt')
		outfile1 = path.join(self.outdir, 'survival1.txt')
		yield 'combine', infile1, False, True, True, dict(arrange = dict(nrow=2, ncol=2)), None, 3, outfile1
		yield 'nc', infile1, False, False, True, dict(params = dict(pval = 'p = {pval}')), None, 3, outfile1
		infile3  = path.join(self.indir, 'survival2.txt')
		cvfile3  = path.join(self.indir, 'survival.cv.txt')
		outfile3 = path.join(self.outdir, 'survival2.txt')
		yield 'cv', infile3, True, False, True, {}, cvfile3, 1, outfile3
		yield 'nop', infile3, True, False, True, False, cvfile3, 1, outfile3

	def testSurvival(self, tag, infile, rnames, combine, pval, plot, covfile, nthread, outfile, devpars = dict(res = 100, height = 300, width = 300)):
		pSurvivalTest                    = pSurvival.copy(tag = tag)
		pSurvivalTest.input              = [infile]
		pSurvivalTest.args.covfile       = covfile
		pSurvivalTest.args.nthread       = nthread
		pSurvivalTest.args.inopts.rnames = rnames
		pSurvivalTest.args.combine       = combine
		#pSurvivalTest.args.devpars       = devpars
		pSurvivalTest.args.pval          = pval
		if isinstance(plot, dict):
			pSurvivalTest.args.plot.update(plot)
		else:
			pSurvivalTest.args.plot = plot
		PyPPL(config).start(pSurvivalTest).run()
		self.assertFileEqual(pSurvivalTest.channel.get(0), outfile)
		
	def dataProvider_testChiSquare(self):
		infile   = path.join(self.indir, 'chisquare1.conttable.txt')
		outfile  = path.join(self.outdir, 'chisquare1.chi2.txt')
		obsvfile = path.join(self.outdir, 'chisquare1.obsv.txt')
		exptfile = path.join(self.outdir, 'chisquare1.expt.txt')
		args     = {'intype': 'cont'}
		yield 't1', infile, outfile, obsvfile, exptfile, args
		yield 't2', path.join(self.indir, 'chisquare2.raw.txt'), outfile, obsvfile, exptfile, {'intype': 'raw'}
		
	def testChiSquare(self, tag, infile, outfile, obsvfile, exptfile, args):
		pChiSquareTest = pChiSquare.copy(tag = tag)
		pChiSquareTest.input = [infile]
		pChiSquareTest.args.update(args)
		PyPPL(config).start(pChiSquareTest).run()
		self.assertFileEqual(pChiSquareTest.channel.get(0), outfile)
		self.assertFileEqual(pChiSquareTest.channel.get(1), obsvfile)
		self.assertFileEqual(pChiSquareTest.channel.get(2), exptfile)
			
	def dataProvider_testFisherExact(self):
		infile   = path.join(self.indir, 'chisquare1.conttable.txt')
		outfile  = path.join(self.outdir, 'fisherexact1.txt')
		args     = {'intype': 'cont'}
		yield 't1', infile, outfile, args
		yield 't2', path.join(self.indir, 'chisquare2.raw.txt'), outfile, {'intype': 'raw'}
		
	def testFisherExact(self, tag, infile, outfile, args):
		pFisherExactTest = pFisherExact.copy(tag = tag)
		pFisherExactTest.input = [infile]
		pFisherExactTest.args.update(args)
		PyPPL(config).start(pFisherExactTest).run()
		self.assertFileEqual(pFisherExactTest.channel.get(0), outfile)
		
	def dataProvider_testPWFisherExact(self):
		yield (
			't1', 
			path.join(self.indir, 'PWFisherExact1.txt'), 
			path.join(self.outdir, 'PWFisherExact1-out.txt'),
			{'intype': 'raw'}
		)
		yield (
			't2', 
			path.join(self.indir, 'PWFisherExact2-pair.txt'), 
			path.join(self.outdir, 'PWFisherExact2-out.txt'),
			{'intype': 'pair'}
		)
		
	def testPWFisherExact(self, tag, infile, outfile, args):
		pPWFisherExactTest = pPWFisherExact.copy(tag = tag)
		pPWFisherExactTest.input = [infile]
		pPWFisherExactTest.args.update(args)
		PyPPL(config).start(pPWFisherExactTest).run()
		self.assertFileEqual(pPWFisherExactTest.channel.get(0), outfile)
	
	def dataProvider_testMediation(self):
		infile = path.join(self.indir, 'mediation.txt')
		exptfile = path.join(self.outdir, 'mediation.txt')
		args = {
			'inopts': {'rnames': False}
		}
		yield infile, args, exptfile
	
	def testMediation(self, infile, args, exptfile):
		pMediation.input = [infile]
		pMediation.args.update(args)
		PyPPL(config).start(pMediation).run()
		self.assertFileEqual(pMediation.channel.get(0), exptfile)

	def dataProvider_testHypergeom(self):
		infile1 = path.join(self.indir, 'hyperg-raw.txt')
		infile2 = path.join(self.indir, 'hyperg-numbers-wN.txt')
		infile3 = path.join(self.indir, 'hyperg-numbers-woN.txt')
		outfile = path.join(self.outdir, 'hyperg.out')
		yield 'raw', infile1, 'raw', 50, outfile
		yield 'wN', infile2, 'numbers', None, outfile
		yield 'woN', infile3, 'numbers', 50, outfile
	
	def testHypergeom(self, tag, infile, intype, N, outfile, inopts = None):
		inopts = inopts or {}
		pHypergeomTest = pHypergeom.copy(tag = tag)
		pHypergeomTest.input = [infile]
		pHypergeomTest.args.intype = intype
		pHypergeomTest.args.inopts = inopts
		pHypergeomTest.args.N = N
		PyPPL(config).start(pHypergeomTest).run()
		self.assertFileEqual(pHypergeomTest.channel.get(0), outfile)
		

if __name__ == '__main__':
	testly.main(failfast = True)