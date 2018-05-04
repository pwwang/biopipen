import unittest, tempfile, shutil, re, mimetypes, sys, testly
from glob import glob
from os import makedirs
from hashlib import md5
from pyppl import Box
from six import StringIO, with_metaclass, assertRaisesRegex as sixAssertRaisesRegex
import inspect, gzip
from subprocess import Popen, PIPE
from pyppl import PyPPL, Proc
from os import path, listdir
from contextlib import contextmanager

filedir = None
config  = {'log': {'file': None, 'level': 'basic', 'lvldiff': ["+P.ARGS", "-DEBUG"]}}

def _getFiledir():
	global filedir
	fn = inspect.getframeinfo(inspect.getouterframes(inspect.currentframe())[2][0])[0]
	fn = path.basename(fn)
	fn = fn[4].lower() + fn[5:-3]
	filedir = path.join(path.realpath(path.dirname(path.dirname(__file__))), 'testfiles', fn)

_getFiledir()

def getfile (name = '', input = True):
	return path.join(filedir, 'input' if input else 'expect', name)

def getbin (name = ''):
	return path.join(path.dirname(path.dirname(path.dirname(filedir))), 'bin', name)

def fileOK(predfile, exptfile, test, order=True, comment=None):
	openfunc = gzip.open if exptfile.endswith('.gz') else open
	with openfunc(predfile) as fpred, openfunc(exptfile) as fexpt:
		predlines = [line for line in fpred.read().splitlines() if not comment or not line.startswith(str(comment))]
		exptlines = [line for line in fexpt.read().splitlines() if not comment or not line.startswith(str(comment))]
		if order:
			test.assertEqual(predlines, exptlines)
		else:
			test.assertEqual(sorted(predlines), sorted(exptlines))

def fileOKIn(exptfile, msg, test):
	openfunc = gzip.open if exptfile.endswith('.gz') else open
	if not isinstance(msg, list):
		msg = [msg]
	with openfunc(exptfile) as f:
		content = f.read()
		for m in msg:
			test.assertIn(m, content)

def procOK(proc, name, test, order = True, comment = None):
	predfile = proc.channel.get()
	exptfile = getfile(name, False)
	if path.isfile(predfile):
		fileOK(predfile, exptfile, test, order, comment)
	else:
		for item in listdir(predfile):
			pfile = path.join(predfile, item)
			efile = path.join(exptfile, item)
			if not path.isfile(pfile): continue
			fileOK(pfile, efile, test, order, comment)

def procOKIn(proc, msg, test):
	fileOKIn(proc.channel.get(), msg, test)

def cmdOK(cmd, test, inout = None, inerr = None, testrc = True):
	p = Popen(cmd, stdout = PIPE, stderr = PIPE)
	out, err = p.communicate()
	if inout:
		if not isinstance(inout, list):
			inout = [inout]
		for ino in inout:
			test.assertIn(ino, out)
	if inerr:
		if not isinstance(inerr, list):
			inerr = [inerr]
		for ine in inerr:
			test.assertIn(ine, err)
	if testrc:
		rc = p.returncode
		test.assertEqual(rc, 0)

def utilTest(input, script, name, tplenvs, test, args = None):
	ends = {
		'.r' : 'Rscript',
		'.py': 'python'
	}
	pTest         = Proc(desc = 'Test utils.', tag=name.split('.')[0])
	pTest.input   = input
	pTest.output  = "outfile:file:outfile"
	pTest.lang    = [ends[k] for k in ends if script.endswith(k)][0]
	pTest.tplenvs = tplenvs
	pTest.args    = {} if not args else args
	pTest.script  = 'file:' + script
	PyPPL(config).start(pTest).run()
	procOK(pTest, name, test)

######

@contextmanager
def captured_output():
	new_out, new_err = StringIO(), StringIO()
	old_out, old_err = sys.stdout, sys.stderr
	try:
		sys.stdout, sys.stderr = new_out, new_err
		yield sys.stdout, sys.stderr
	finally:
		sys.stdout, sys.stderr = old_out, old_err

def md5sum(fn):
	ret = md5()
	with open(fn, "rb") as f:
		#for chunk in iter(lambda: f.read(4096), b""):
		#	ret.update(chunk)
		ret.update(f.read())
	return ret.hexdigest()

class TestCase(testly.TestCase):

	def assertItemEqual(self, first, second, msg = None):
		first          = [repr(x) for x in first]
		second         = [repr(x) for x in second]
		first          = str(sorted(first))
		second         = str(sorted(second))
		assertion_func = self._getAssertEqualityFunc(first, second)
		assertion_func(first, second, msg=msg)

	def assertDictIn(self, first, second, msg = 'Not all k-v pairs in 1st element are in the second.'):
		assert isinstance(first, dict)
		assert isinstance(second, dict)
		notInkeys = [k for k in first.keys() if k not in second.keys()]
		if notInkeys:
			self.fail(msg = 'Keys of first dict not in second: %s' % notInkeys)
		else:
			seconds2 = {k:second[k] for k in first.keys()}
			for k in first.keys():
				v1   = first[k]
				v2   = second[k]
				try:
					self.assertSequenceEqual(v1, v2)
				except AssertionError:
					self.assertEqual(v1, v2)



	def assertDictNotIn(self, first, second, msg = 'all k-v pairs in 1st element are in the second.'):
		assert isinstance(first, dict)
		assert isinstance(second, dict)
		ret = False
		for k in first.keys():
			if k in second:
				if first[k] != second[k]:
					ret = True
			else:
				ret = True
		if not ret:
			self.fail(msg)

	def assertFileEqual(self, first, second, msg = None):
		import icdiff
		import filecmp
		import sys
		if not filecmp.cmp(first, second, shallow=False):
			maxdiff = self.maxDiff or 0
			origargv = [arg for arg in sys.argv]
			sys.argv[1:] = ['-L', first, '-L', second, '-U', '1', '--head', str(int(maxdiff))]
			options = icdiff.get_options()[0]
			sys.argv = origargv
			icdiff.diff(first, second, options)
			msg = msg or ''
			self.fail('%s\nSet self.maxDiff = None to see all diff.' % msg)

	def assertFileCountEqual(self, first, second, sort = 'sort', msg = None):
		tmpdir = tempfile.gettempdir()
		firstsorted  = path.join(tmpdir, path.splitext(path.basename(first))[0] + '.sorted')
		secondsorted = path.join(tmpdir, path.splitext(path.basename(second))[0] + '.sorted')
		Popen([sort, '-o', firstsorted, first]).wait()
		Popen([sort, '-o', secondsorted, second]).wait()
		import icdiff
		import filecmp
		import sys
		if not filecmp.cmp(firstsorted, secondsorted, shallow=False):
			maxdiff = self.maxDiff or 0
			origargv = [arg for arg in sys.argv]
			sys.argv[1:] = ['-L', first, '-L', second, '-U', '1', '--head', str(int(maxdiff))]
			options = icdiff.get_options()[0]
			sys.argv = origargv
			icdiff.diff(firstsorted, secondsorted, options)
			msg = msg or ''
			self.fail('%s\nSet self.maxDiff = None to see all diff.' % msg)


	def assertDirEqual(self, first, second, msg = None):
		if not path.isdir(first):
			self.fail('The first file is not a directory.')
		if not path.isdir(second):
			self.fail('The second file is not a directory.')
		for fn in glob(path.join(first, '*')):
			bn = path.basename(fn)
			if path.isdir(fn):
				self.assertDirEqual(fn, path.join(second, bn))
			else:
				self.assertFileEqual(fn, path.join(second, bn))

	def assertTextEqual(self, first, second, msg = None):
		if not isinstance(first, list):
			first  = first.split('\n')
		if not isinstance(second, list):
			second = second.split('\n')
		self.assertListEqual(first, second, msg)

	def assertRaisesStr(self, exc, s, callable, *args, **kwds):
		sixAssertRaisesRegex(self, exc, s, callable, *args, **kwds)

	def assertItemSubset(self, s, t, msg = 'The first list is not a subset of the second.'):
		assert isinstance(s, list)
		assert isinstance(t, list)
		self.assertTrue(set(s) < set(t), msg = msg)

	def assertInFile(self, s, f):
		sf = readFile(f, str)
		self.assertIn(s, sf)

	def assertBamEqual(self, bam1, bam2, samtools = 'samtools', msg = None):
		tmpdir = tempfile.gettempdir()
		sam1   = path.join(tmpdir, path.splitext(bam1)[0] + '.sam')
		sam2   = path.join(tmpdir, path.splitext(bam2)[0] + '.sam')
		Popen([samtools, 'view', '-h', '-o', sam1, bam1]).wait()
		Popen([samtools, 'view', '-h', '-o', sam2, bam2]).wait()
		import icdiff
		import filecmp
		import sys
		if not filecmp.cmp(sam1, sam2, shallow=False):
			maxdiff = self.maxDiff or 0
			origargv = [arg for arg in sys.argv]
			sys.argv[1:] = ['-L', bam1, '-L', bam2, '-U', '1', '--head', str(int(self.maxDiff))]
			options = icdiff.get_options()[0]
			sys.argv = origargv
			icdiff.diff(sam1, sam2, options)
			msg = msg or ''
			self.fail('%s\nSet self.maxDiff = None to see all diff.' % msg)

	def assertBamCountEqual(self, bam1, bam2, samtools = 'samtools', msg = None):
		tmpdir = tempfile.gettempdir()
		sam1   = path.join(tmpdir, path.splitext(bam1)[0] + '.sam')
		sam2   = path.join(tmpdir, path.splitext(bam2)[0] + '.sam')
		Popen([samtools, 'sort', '-n', '-o', sam1, bam1]).wait()
		Popen([samtools, 'sort', '-n', '-o', sam2, bam2]).wait()
		import icdiff
		import filecmp
		import sys
		if not filecmp.cmp(sam1, sam2, shallow=False):
			maxdiff = self.maxDiff or 0
			origargv = [arg for arg in sys.argv]
			sys.argv[1:] = ['-L', bam1, '-L', bam2, '-U', '1', '--head', str(int(self.maxDiff))]
			options = icdiff.get_options()[0]
			sys.argv = origargv
			icdiff.diff(sam1, sam2, options)
			msg = msg or ''
			self.fail('%s\nSet self.maxDiff = None to see all diff.' % msg)

def testdirs(classname):
	testdir   = path.join(tempfile.gettempdir(), 'bioprocs_unittest', classname)
	if path.exists(testdir):
		shutil.rmtree(testdir)
	makedirs(testdir)
	parentdir = path.join(path.dirname(path.dirname(__file__)), 'testfiles', classname[4].lower() + classname[5:])
	indir     = path.join(parentdir, 'input')
	outdir    = path.join(parentdir, 'expect')
	return testdir, indir, outdir
