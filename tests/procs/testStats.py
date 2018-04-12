import helpers, unittest
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.stats import pMetaPval, pMetaPval1, pSurvival, pChiSquare, pFisherExact, pPWFisherExact, pMediation
from bioprocs.common import pFiles2Dir
from bioprocs.tsv import pRbind

class testStats (helpers.TestCase):
	
	def dataProvider_test1MetaPval(self, testdir, indir, outdir):
		yield 't1', indir, path.join(outdir, 'metapval.txt'), {'pattern': "testMetaPval*.txt"}

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

	def testSurvival(self):
		pSurvival1 = pSurvival.copy()
		pSurvival1.input                = [getfile('survival.txt')]
		pSurvival1.args.rnames          = False
		pSurvival1.args.nthread         = 3
		pSurvival1.args.gridParams.nrow = 2
		pSurvival1.args.gridParams.ncol = 2
		PyPPL(config).start(pSurvival1).run()
		procOK(pSurvival1, 'survival-combine', self)

	def dataProvider_testSurvivalNoCombine(self, testdir, indir, outdir):
		infile  = path.join(indir, 'survival.txt')
		outfile = path.join(outdir, 'survival-individual')
		yield 't1', infile, outfile
	
	def testSurvivalNoCombine(self, tag, infile, outfile):
		pSurvival2              = pSurvival.copy(tag = tag)
		pSurvival2.input        = [infile]
		pSurvival2.args.rnames  = False
		pSurvival2.args.nthread = 1
		pSurvival2.args.combine = False
		PyPPL(config).start(pSurvival2).run()
		self.assertDirEqual(pSurvival2.channel.get(), outfile)
		
	def dataProvider_testChiSquare(self, testdir, indir, outdir):
		infile   = path.join(indir, 'chisquare1.conttable.txt')
		outfile  = path.join(outdir, 'chisquare1.chi2.txt')
		obsvfile = path.join(outdir, 'chisquare1.obsv.txt')
		exptfile = path.join(outdir, 'chisquare1.expt.txt')
		args     = {'intype': 'cont'}
		yield 't1', infile, outfile, obsvfile, exptfile, args
		yield 't2', path.join(indir, 'chisquare2.raw.txt'), outfile, obsvfile, exptfile, {'intype': 'raw'}
		
	def testChiSquare(self, tag, infile, outfile, obsvfile, exptfile, args):
		pChiSquareTest = pChiSquare.copy(tag = tag)
		pChiSquareTest.input = [infile]
		pChiSquareTest.args.update(args)
		PyPPL(config).start(pChiSquareTest).run()
		self.assertFileEqual(pChiSquareTest.channel.get(0), outfile)
		self.assertFileEqual(pChiSquareTest.channel.get(1), obsvfile)
		self.assertFileEqual(pChiSquareTest.channel.get(2), exptfile)
			
	def dataProvider_testFisherExact(self, testdir, indir, outdir):
		infile   = path.join(indir, 'chisquare1.conttable.txt')
		outfile  = path.join(outdir, 'fisherexact1.txt')
		args     = {'intype': 'cont'}
		yield 't1', infile, outfile, args
		yield 't2', path.join(indir, 'chisquare2.raw.txt'), outfile, {'intype': 'raw'}
		
	def testFisherExact(self, tag, infile, outfile, args):
		pFisherExactTest = pFisherExact.copy(tag = tag)
		pFisherExactTest.input = [infile]
		pFisherExactTest.args.update(args)
		PyPPL(config).start(pFisherExactTest).run()
		self.assertFileEqual(pFisherExactTest.channel.get(0), outfile)
		
	def dataProvider_testPWFisherExact(self, testdir, indir, outdir):
		yield (
			't1', 
			path.join(indir, 'PWFisherExact1.txt'), 
			path.join(outdir, 'PWFisherExact1-out.txt'),
			{'intype': 'raw'}
		)
		yield (
			't2', 
			path.join(indir, 'PWFisherExact2-pair.txt'), 
			path.join(outdir, 'PWFisherExact2-out.txt'),
			{'intype': 'pair'}
		)
		
	def testPWFisherExact(self, tag, infile, outfile, args):
		pPWFisherExactTest = pPWFisherExact.copy(tag = tag)
		pPWFisherExactTest.input = [infile]
		pPWFisherExactTest.args.update(args)
		PyPPL(config).start(pPWFisherExactTest).run()
		self.assertFileEqual(pPWFisherExactTest.channel.get(0), outfile)
	
	def dataProvider_testMediation(self, testdir, indir, outdir):
		infile = path.join(indir, 'mediation.txt')
		exptfile = path.join(outdir, 'mediation.txt')
		args = {
			'inopts': {'rnames': False}
		}
		yield infile, args, exptfile
	
	def testMediation(self, infile, args, exptfile):
		pMediation.input = [infile]
		pMediation.args.update(args)
		PyPPL(config).start(pMediation).run()
		self.assertFileEqual(pMediation.channel.get(0), exptfile)
		

if __name__ == '__main__':
	unittest.main(failfast = True)