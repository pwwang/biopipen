import testly
from helpers import runBioprocs, getData, getOutput


class TestAlgorithm(testly.TestCase):

	def dataProvider_testAR(self):
		yield getData('algorithm.pAR.paper', **{
			'args.seed' : 8525,
			'args.tfrac': .472,
			'args.svdP' : 25 # may have some issues with svd in R
		}) , getData('algorithm.pAR.paper', 'e')
		yield getData('algorithm.pAR.default'), getData('algorithm.pAR.default', 'e')

	def testAR(self, args, expectfiles = None):
		c = runBioprocs('algorithm.pAR', args)
		expectfiles = expectfiles or {}
		outfiles = getOutput(c)
		for key, val in expectfiles.items():
			k = key.split('.', 1)[-1]
			self.assertFileEqual(val, outfiles[k])
		

if __name__ == '__main__':
	testly.main(verbosity = 2)

