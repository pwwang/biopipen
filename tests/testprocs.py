import testly
from collections import defaultdict
from helpers import runBioprocs, getData, getOutput

def TestClassFactory(klass, datarows):

	def dataprovider_maker(self):
		for datarow in datarows:
			yield datarow

	def testit(self, proc, args, exptfiles = None, opt1 = None, opt2 = None):
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

	def __init__(self, *args, **kwargs):
		testly.TestCase.__init__(self, *args, **kwargs)

	members = {'__init__': __init__}
	proc = datarows[0].kwargs['proc'][len(klass)+1:]
	proc = 'test' + proc
	dataprovider = 'dataProvider_' + proc
	members[dataprovider] = dataprovider_maker
	members[proc] = testit

	globals()['Test' + klass[0].upper() + klass[1:]] = type(klass, (testly.TestCase,), members)

DATAROWS = defaultdict(lambda: [])
for data in getData():
	klass = data.kwargs['proc'].split('.')[0]
	DATAROWS[klass].append(data)

for key, val in DATAROWS.items():
	TestClassFactory(key, val)

if __name__ == '__main__':
	testly.main(verbosity = 2)