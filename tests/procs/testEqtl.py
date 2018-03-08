import helpers, unittest

from os import path
from pyppl import PyPPL
from bioprocs.eqtl import pMatrixeQTL

class TestEqtl(helpers.TestCase):
	
	def dataProvider_testMatrixeQTL(self, testdir, indir, outdir):
		snpfile    = path.join(indir, 'SNP.txt')
		snpfile_na = path.join(indir, 'SNP_NA.txt')
		expfile    = path.join(indir, 'GE.txt')
		notfound   = path.join(outdir, 'notfound.txt')
		e_2        = path.join(outdir, '1e-2.txt')
		allparams  = path.join(outdir, 'allparams.txt')
		snpna      = path.join(outdir, 'snpNA.txt')
		snppos     = path.join(indir, 'snpsloc.txt')
		genepos    = path.join(indir, 'geneloc.txt')
		covfile    = path.join(indir, 'Covariates.txt')
		yield snpfile, expfile, notfound
		yield snpfile, expfile, e_2, '', {'pval': 1e-2}
		yield snpfile_na, expfile, snpna, '', {'pval': 2e-2}
		yield snpfile, expfile, allparams, covfile, {'pval': 2e-2, 'cisopts': {'snppos': snppos, 'genepos': genepos, 'dist': 1e6, 'cispv': 3e-2}}
		
	
	def testMatrixeQTL(self, snpfile, expfile, outfile, covfile = '', args = None):
		args = args or {}
		pMatrixeQTL1 = pMatrixeQTL.copy(tag = 'tag' + self.index())
		pMatrixeQTL1.args.update(args)
		pMatrixeQTL1.input = snpfile, expfile, covfile
		PyPPL(helpers.config).start(pMatrixeQTL1).run()
		self.assertFileEqual(pMatrixeQTL1.channel.get(), outfile)
	
	
if __name__ == '__main__':
	unittest.main()
		

