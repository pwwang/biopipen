import helpers, unittest, testly

from os import path
from pyppl import PyPPL
from bioprocs.eqtl import pMatrixeQTL

class TestEqtl(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestEqtl')

	def dataProvider_testMatrixeQTL(self):
		snpfile    = path.join(self.indir, 'SNP.txt')
		snpfile_na = path.join(self.indir, 'SNP_NA.txt')
		expfile    = path.join(self.indir, 'GE.txt')
		notfound   = path.join(self.outdir, 'notfound.txt')
		e_2        = path.join(self.outdir, '1e-2.txt')
		allparams  = path.join(self.outdir, 'allparams.txt')
		snpna      = path.join(self.outdir, 'snpNA.txt')
		snppos     = path.join(self.indir, 'snpsloc.txt')
		genepos    = path.join(self.indir, 'geneloc.txt')
		covfile    = path.join(self.indir, 'Covariates.txt')
		yield 't1', snpfile, expfile, notfound
		yield 't2', snpfile, expfile, e_2, '', {'pval': 1e-2}
		yield 't3', snpfile_na, expfile, snpna, '', {'pval': 2e-2}
		yield 't4', snpfile, expfile, allparams, covfile, {'pval': 2e-2, 'cisopts': {'snppos': snppos, 'genepos': genepos, 'dist': 1e6, 'cispv': 3e-2}}


	def testMatrixeQTL(self, tag, snpfile, expfile, outfile, covfile = '', args = None):
		args = args or {}
		pMatrixeQTL1 = pMatrixeQTL.copy(tag = tag)
		pMatrixeQTL1.args.update(args)
		pMatrixeQTL1.input = snpfile, expfile, covfile
		PyPPL(helpers.config).start(pMatrixeQTL1).run()
		self.assertFileEqual(pMatrixeQTL1.channel.get(), outfile)


if __name__ == '__main__':
	testly.main()
