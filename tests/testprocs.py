import testly
from helpers import runBioprocs, getData, getOutput

class TestProcs(testly.TestCase):

	def dataProvider_testProc(self):
		yield 'mlearn.pCrossValid', getData('mlearn.pCrossValid.default', **{
			'args.seed'        : 998,
			'args.train.method': 'rf',
			'args.train.form'  : 'Class ~ .',
			'args.ctrl.method' : 'cv',
			'args.ctrl.number' : 10
		}) , getData('mlearn.pCrossValid.default', 'e')

		yield 'algorithm.pAR', getData('algorithm.pAR.paper', **{
			'args.seed' : 8525,
			'args.tfrac': .472,
			'args.svdP' : 25 # may have some issues with svd in R
		}) , getData('algorithm.pAR.paper', 'e')

		yield 'algorithm.pAR', getData('algorithm.pAR.default', **{
			'args.seed': 8525
		}), getData('algorithm.pAR.default', 'e')

		yield 'snp.pRs2Bed', getData('snp.pRs2Bed.default', **{
			'args.inopts.cnames': False
		}) , getData('snp.pRs2Bed.default', 'e')



	def testProc(self, proc, args, exptfiles = None):
		c = runBioprocs(proc, args)
		exptfiles = exptfiles or {}
		outfiles = getOutput(c)
		for key, val in exptfiles.items():
			k = key.split('.', 1)[-1]
			if k in outfiles:
				self.assertFileEqual(val, outfiles[k])
			else:
				self.assertFileEqual(val, k.format(**outfiles))

if __name__ == '__main__':
	testly.main(verbosity = 2)