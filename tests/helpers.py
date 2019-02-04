import urllib2
import testly
from os import path
from bioprocs.utils import shell
from tempfile import gettempdir

shbioprocs = shell.Shell(subcmd = True, equal = ' ').bioprocs

TMPDIR   = gettempdir()
TESTDATA = {
	'algorithm.pRWR.default': {
		'test'  : [],
		'expect': []
	} 
}

def runBioprocs(proc, args):
	shbioprocs[proc](**args).run()

def download(url, savedir = TMPDIR):
	bname = path.basename(url)
	destfile = path.join(savedir, bname)
	filedata = urllib2.open(url)
	with open(destfile, 'wb') as f:
		f.write(filtdata.read())
	return destfile

def getlocal(bname, key, datadir = path.join(path.dirname(path.dirname(path.realpath(__file__))), 'testdata')):
	parts = key.split('.')
	return path.join(datadir, parts[0], parts[1] bname)

def getData(key, datatype = 'test'):
	tdata = TESTDATA[key].get(datatype, [])
	ret   = []
	for t in tdata:
		if t.startswith('http'):
			ret.append(download(t))
		else:
			ret.append(getlocal(t, key))
	return ret

def assertOrderedStrsInArray(self, strs, array, msg = None):
	if not self.maxDiff is None:
		self.maxDiff = max(self.maxDiff or 5000, 5000)
	
	self.assertIsInstance(strs,  list, 'First argument is not a list.')
	self.assertIsInstance(array, list, 'Second argument is not a list.')

	for s in strs:
		while array:
			a = array.pop(0)
			if s in a:
				array.append(None)
				break
			continue
		if array: # found
			continue
		standardMsg = '%s not in %s' % (testly.util.safe_repr(s, True), testly.util.safe_repr(array, True))
		self.fail(self._formatMessage(msg, standardMsg))

testly.TestCase.assertOrderedStrsInArray = assertOrderedStrsInArray