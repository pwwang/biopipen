import re
import urllib2
import testly
import yaml
from os import path
from bioprocs.utils import shell
from bioprocs.utils.tsvio2 import TsvReader
from tempfile import gettempdir

shbioprocs = shell.Shell(subcmd = True, dash = '-', equal = ' ', duplistkey = False).bioprocs

PROCDATADIR = path.join(path.dirname(path.abspath(__file__)), 'procdata')
DATAFILE    = path.join(PROCDATADIR, 'data.yml')
TMPDIR      = gettempdir()
CACHED      = {}

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
	d1, d2 = key.split('.')
	return path.join(datadir, d1, d2, bname)

def getioval(value, proc = None):
	if isinstance(value, list):
		return [getioval(v, proc) for v in value]
	if not value.startswith('plain:') and not value.startswith('http') and \
	not value.startswith('ftp') and not value.startswith('file:'):
		value = 'file:' + value
	if value in CACHED:
		return CACHED[value]
	elif value.startswith('plain:'):
		return value[7:]
	elif value.startswith('http') or value.startswith('ftp'):
		CACHED[value] = download(v)
		return CACHED[value]
	elif value.startswith('file:'):
		return getlocal(value[5:], proc)
	return value

def getData(datafile = DATAFILE):
	with open(datafile) as stream:
		data = yaml.load(stream)

	for key, value in data.items():
		kparts    = key.split('.')
		tag       = kparts[-1] if len(kparts) > 2 else 'default'
		proc      = '.'.join(kparts[:2])
		args      = {'tag': tag}
		exptfiles = {}
		opt1      = {}
		opt2      = {}
		for k, v in value.items():
			ret = args
			if k[:2] in ['i.', 'o.']:
				v = getioval(v, proc)
				if k[:2] == 'o.':
					k = k[2:]
					exptfiles[k] = v
				elif isinstance(v, list):
					args[k + ':list:one'] = v
				else:
					args[k] = v
			elif k.startswith('opt'):
				ret = opt1 if k.startswith('opt1.') else opt2
				ret[k[5:]] = v
			else:
				ret[k] = v
		yield testly.Data(proc = proc, args = args, exptfiles = exptfiles, opt1 = opt1, opt2 = opt2)

def getOutput(c):
	ret = {}
	for line in c.stderr.splitlines():
		#[2019-02-06 09:49:48  OUTPUT] pAR: [1/1] outdir => /local2/tmp/m161047/bioprocs.workdir/PyPPL.pAR.notag.6duu519e/1/output/motif_hits-protein_expression_t-gene_expression.AR
		m = re.match(r'^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}  OUTPUT\] \w+: \[\d+\/\d+\] (\w+) \s*=> (.+)$', line)
		if not m: continue
		ret[m.group(1)] = m.group(2)
	return ret


def istext(filename):
	import string 
	with open(filename) as f:
		s = f.read(512)
	text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
	_null_trans = string.maketrans("", "")
	if not s:
		# Empty files are considered text
		return True
	if "\0" in s:
		# Files with null bytes are likely binary
		return False
	# Get the non-text characters (maps a character to itself then
	# use the 'remove' option to get rid of the text characters.)
	t = s.translate(_null_trans, text_characters)
	# If more than 30% non-text characters, then
	# this is considered a binary file
	if float(len(t))/float(len(s)) > 0.30:
		return False
	return True


# unittest asserts
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

def assertFileEqual(self, first, second, filetype = None, firstInopts = None, secondInopts = None, msg = None):
	if not self.maxDiff is None:
		self.maxDiff = max(self.maxDiff or 5000, 5000)

	filetype1 = filetype or ('text' if istext(first) else 'nontext')
	filetype2 = filetype or ('text' if istext(second) else 'nontext')
	if filetype1 != filetype2:
		standardMsg = 'Files different, because file1 is {0} but file2 is {1}'.format(
			filetype1, filetype2
		)
		self.fail(self._formatMessage(msg, standardMsg))
	elif filetype1 == 'text':# and filetype2 == 'text':
		reader1 = TsvReader(first,  **firstInopts)  if firstInopts  else TsvReader(first)
		reader2 = TsvReader(second, **secondInopts) if secondInopts else TsvReader(second)
		rindex  = 0
		for r1 in reader1:
			rindex += 1
			try:
				r2 = next(reader2)
			except StopIteration:
				standardMsg = 'File1 and file2 are different.\nFile1: {2}\nFile2: {3}\nRow {0} of file1 is: {1}, but nothing at row {0} of file2.'.format(rindex, r1, first, second)
				self.fail(self._formatMessage(msg, standardMsg))
			if r1 != r2:
				standardMsg = 'File1 and file2 are different.\nFile1: {3}\nFile2: {4}\nRow {0} of file1: {1}\nRow {0} of file2: {2}'.format(rindex, r1, r2, first, second)
				self.fail(self._formatMessage(msg, standardMsg))
	else: # filetype1 == 'nontext' and filetype2 == 'nonetext': # binary
		import filecmp
		if not filecmp.cmp(first, second, shallow = False):
			standardMsg = 'Binary files are different:\n{}\n{}'.format(first, second)
			self.fail(self._formatMessage(msg, standardMsg))

testly.TestCase.assertOrderedStrsInArray = assertOrderedStrsInArray
testly.TestCase.assertFileEqual = assertFileEqual