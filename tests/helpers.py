import re
import urllib2
import testly
from os import path
from bioprocs.utils import shell, alwaysList
from bioprocs.utils.tsvio2 import TsvReader
from tempfile import gettempdir

shbioprocs = shell.Shell(subcmd = True, dash = '-', equal = ' ').bioprocs

PROCDATADIR = path.join(path.dirname(path.dirname(path.abspath(__file__))), 'procdata')
TMPDIR      = gettempdir()
TESTDATA    = {
	'algorithm.pRWR.default': {
		'tag'   : 'default',
		'test'  : [],
		'expect': []
	},
	'algorithm.pAR.default': {
		'tag' : 'default',
		'i.D' : 'ar.d.txt.gz',
		'i.Pt': 'ar.pt.txt.gz',
		'i.Y' : 'ar.y.txt.gz',
		'e.W' : 'default.W.txt',
	},
	'algorithm.pAR.paper': {
		'tag' : 'paper',
		'i.D' : 'motif_hits.txt',
		'i.Pt': 'protein_expression_t.txt',
		'i.Y' : 'gene_expression.txt',
		'e.W' : 'paper.W.txt',
	} 
}

def runBioprocs(proc, args):
	args['config._log.shortpath:py'] = 'False'
	return shbioprocs[proc](**args).run(save = 'same', uselogger = False)

def download(url, savedir = TMPDIR):
	bname = path.basename(url)
	destfile = path.join(savedir, bname)
	filedata = urllib2.open(url)
	with open(destfile, 'wb') as f:
		f.write(filtdata.read())
	return destfile

def getlocal(bname, key, datadir = PROCDATADIR):
	parts = key.split('.')
	return path.join(datadir, parts[0], parts[1], bname)

def getData(key, which = 'i., tag', **kwargs):
	which = alwaysList(which)
	ret   = kwargs
	tdata = TESTDATA[key]
	for k, v in tdata.items():
		if not any(k.startswith(w) for w in which): continue
					
		if v.startswith('http'):
			data = download(v)
			TESTDATA[key][k] = 'cached:' + data
			ret[k] = data
		elif v.startswith('cached:'):
			ret[k] = v[7:]
		elif k.startswith('i.') or k.startswith('e.'):
			ret[k] = getlocal(v, key)
		else:
			ret[k] = v
	return ret

def getOutput(c):
	ret = {}
	for line in c.stderr.splitlines():
		#[2019-02-06 09:49:48  OUTPUT] pAR: [1/1] outdir => /local2/tmp/m161047/bioprocs.workdir/PyPPL.pAR.notag.6duu519e/1/output/motif_hits-protein_expression_t-gene_expression.AR
		m = re.match(r'^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}  OUTPUT\] \w+: \[\d+\/\d+\] (\w+) \s*=> (.+)$', line)
		if not m: continue
		ret[m.group(1)] = m.group(2)
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

def assertFileEqual(self, first, second, firstInopts = None, secondInopts = None, msg = None):
	if not self.maxDiff is None:
		self.maxDiff = max(self.maxDiff or 5000, 5000)

	reader1 = TsvReader(first,  firstInopts)  if firstInopts  else TsvReader(first)
	reader2 = TsvReader(second, secondInopts) if secondInopts else TsvReader(second)
	rindex  = 0
	for r1 in reader1:
		rindex += 1
		try:
			r2 = next(reader2)
		except StopIteration:
			standardMsg = 'Row {0} of file1 is: {1}, but nothing at row {0} of file2.'.format(rindex, r1)
			self.fail(self._formatMessage(msg, standardMsg))
		if r1 != r2:
			standardMsg = 'Row {0} of file1: {1}\nRow {0} of file2: {2}'.format(rindex, r1, r2)
			self.fail(self._formatMessage(msg, standardMsg))

testly.TestCase.assertOrderedStrsInArray = assertOrderedStrsInArray
testly.TestCase.assertFileEqual = assertFileEqual