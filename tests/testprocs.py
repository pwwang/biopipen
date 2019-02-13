import testly
from helpers import runBioprocs, getData, getOutput

class TestProcs(testly.TestCase):

	def dataProvider_testProc(self):
		for data in getData():
			yield data

	def testProc(self, proc, args, exptfiles = None, opt1 = None, opt2 = None):
		c         = runBioprocs(proc, args)
		exptfiles = exptfiles or {}
		outfiles  = getOutput(c)
		opt1      = opt1 or {}
		opt2      = opt2 or {}
		for key, val in exptfiles.items():
			if key in outfiles:
				self.assertFileEqual(val, outfiles[key], firstInopts = opt1, secondInopts = opt2)
			else:
				self.assertFileEqual(val, key.format(**outfiles), firstInopts = opt1, secondInopts = opt2)

if __name__ == '__main__':
	testly.main(verbosity = 2)