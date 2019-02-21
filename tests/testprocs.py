from __future__ import print_function
import sys
import testly
from collections import defaultdict
from helpers import runBioprocs, getData, getOutput

def TestClassFactory(klass, datarows, listtests = False):

	def dataprovider_maker(drows):
		def df(self):
			for drow in drows:
				yield drow
		return df

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
	procs   = defaultdict(lambda: [])
	for datarow in datarows:
		procs[datarow.kwargs['proc'][len(klass)+1:]].append(datarow)

	for proc, drows in procs.items():
		proc = 'test' + proc
		dataprovider = 'dataProvider_' + proc
		members[dataprovider] = dataprovider_maker(drows)
		members[proc] = testit

	if listtests:
		print('Test' + klass[0].upper() + klass[1:])
		for member in members:
			if member.startswith('__') or member.startswith('dataProvider_'):
				continue
			print('- ' + member)
		print('')
	else:
		globals()['Test' + klass[0].upper() + klass[1:]] = type(klass, (testly.TestCase,), members)

DATAROWS = defaultdict(lambda: [])
for data in getData():
	klass = data.kwargs['proc'].split('.')[0]
	DATAROWS[klass].append(data)

if __name__ == '__main__':
	listtests = len(sys.argv) > 1 and sys.argv[1] == 'list'
	for key, val in DATAROWS.items():
		TestClassFactory(key, val, listtests)

	if not listtests:
		testly.main(verbosity = 2)
