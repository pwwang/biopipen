import helpers, unittest, testly
from os import path
from pyppl import PyPPL
from bioprocs.power import pSurvivalPower

class TestPower (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestPower')

	def dataProvider_testSurvivalPower(self):
		yield (
			't1',
			path.join(self.indir, 'survpower1.txt'),
			path.join(self.outdir, 'survpower1.txt'),
			{'intype': 'ratio', 'plot': True}
		)
		yield (
			't2',
			path.join(self.indir, 'survpower2.txt'),
			path.join(self.outdir, 'survpower2.txt'),
			{'intype': 'detailed', 'plot': True}
		)
		yield (
			't3',
			path.join(self.indir, 'survpower3.txt'),
			path.join(self.outdir, 'survpower3.txt'),
			{'intype': 'detailed', 'plot': True}
		)

	def testSurvivalPower(self, tag, infile, outfile, args = None):
		pSurvivalPowerTest = pSurvivalPower.copy(tag = tag)
		pSurvivalPowerTest.args.update(args or {})
		pSurvivalPowerTest.input = [infile]
		PyPPL(helpers.config).start(pSurvivalPowerTest).run()
		self.assertFileEqual(path.join(pSurvivalPowerTest.channel.get(), path.splitext(path.basename(infile))[0] + '.ssizes'), outfile)



if __name__ == '__main__':
	testly.main()
