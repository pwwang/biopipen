import helpers, unittest
from collections import OrderedDict
from os import path
from subprocess import check_output
from bioprocs.utils.helpers import runcmd, cmdargs, RuncmdException, mem2, Mem2Exception


class TestHelpers(helpers.TestCase):
	
	def dataProvider_testRuncmd(self):
		yield 'ls', True, False
		yield 'cmdnotexists', False, True, RuncmdException
		yield 'cmdnotexists', False, False
		
	def testRuncmd(self, cmd, ret, quit, exception = None):
		if exception:
			with helpers.captured_output() as (out, err):
				self.assertRaises(exception, runcmd, cmd, quit)
			errvalue = err.getvalue()
			self.assertIn('Running command at PID:', errvalue)
			self.assertIn(cmd, errvalue)
			self.assertIn('Return code: ', errvalue)
		else:
			with helpers.captured_output() as (out, err):
				r = runcmd(cmd, quit)
			self.assertEqual(r, ret)
			errvalue = err.getvalue()
			self.assertIn('Running command at PID:', errvalue)
			self.assertIn(cmd, errvalue)
			self.assertIn('Return code: ', errvalue)
	
	def dataProvider_testRuncmdR(self, testdir, indir):
		srcfile = path.join(indir, 'runcmd1.r')
		yield srcfile, True
		srcfile = path.join(indir, 'runcmd2.r')
		yield srcfile, False
	
	def testRuncmdR(self, srcfile, ret):
		with helpers.captured_output() as (out, err):
			r = runcmd('Rscript ' + srcfile, quit = False)
		self.assertEqual(r, ret)
		
	def dataProvider_testCmdargs(self):
		yield {}, 'auto', 'auto', ''
		yield OrderedDict([
			('c', True),
			('de', 'de fg'),
			('b', 2),
			('a', 1),
		]), 'auto', 'auto', "-c --de='de fg' -b 2 -a 1"
		yield OrderedDict([
			('c', True),
			('de', 'de fg'),
			('b', 2),
			('a', 1),
		]), '--', 'auto', "--c --de='de fg' --b 2 --a 1"
		yield OrderedDict([
			('c', True),
			('de', 'de fg'),
			('b', 2),
			('a', 1),
		]), '--', ' ', "--c --de 'de fg' --b 2 --a 1"
		
	def testCmdargs(self, params, dash, equal, ret):
		r = cmdargs(params, dash = dash, equal = equal, )
		self.assertEqual(r, ret)
		
	def dataProvider_testCmdargsR(self, testdir, indir):
		yield path.join(indir, 'cmdargs1.r'), ''
		yield path.join(indir, 'cmdargs2.r'), "-a 1.0 -c --de='de fg' -b 2.0"
		yield path.join(indir, 'cmdargs3.r'), "--a 1.0 --c --de='de fg' --b 2.0"
		yield path.join(indir, 'cmdargs4.r'), "--a 1.0 --c --de 'de fg' --b 2.0"
		
	def testCmdargsR(self, rfile, ret):
		out = check_output(['Rscript', rfile])
		self.assertEqual(out, ret)
		
	def dataProvider_testMem2(self):
		yield '1x', 'auto', None, Mem2Exception
		yield '123k', 'auto', '123K'
		yield '100M', 'auto', '100M'
		yield '100M', 'K', '102400K'
		yield '100M', 'Java', '-Xms12800K -Xmx100M'
		yield '8G', 'Java', '-Xms1G -Xmx8G'
		
	def testMem2(self, mem, unit, ret, exception = None):
		if exception:
			self.assertRaises(exception, mem2, mem, unit)
		else:
			m = mem2(mem, unit)
			self.assertEqual(m, ret)
	
	def dataProvider_testMem2R(self, testdir, indir):
		rfile = path.join(indir, 'mem2.r')
		yield rfile, '123k', 'auto', '123K'
		yield rfile, '100M', 'auto', '100M'
		yield rfile, '100M', 'K', '102400K'
		yield rfile, '100M', 'Java', '-Xms12800K -Xmx100M'
		yield rfile, '8G', 'Java', '-Xms1G -Xmx8G'
	
	def testMem2R(self, rfile, mem, unit, ret):
		out = check_output(['Rscript', rfile, mem, unit])
		self.assertEqual(out, ret)
		
	def dataProvider_testCbindfill(self, testdir, indir, outdir):
		rfile    = path.join(indir, 'cbindfill.r')
		infile1  = path.join(indir, 'cbindfill1_1.txt')
		infile2  = path.join(indir, 'cbindfill1_2.txt')
		outfile  = path.join(testdir, 'cbindfill1.txt')
		exptfile = path.join(outdir, 'cbindfill1.txt')
		yield rfile, infile1, infile2, outfile, exptfile
		
	def testCbindfill(self, rfile, infile1, infile2, outfile, exptfile):
		with helpers.captured_output():
			runcmd(['Rscript', rfile, infile1, infile2, outfile])
		self.assertFileEqual(outfile, exptfile)
	
	def dataProvider_testRbindfill(self, testdir, indir, outdir):
		rfile    = path.join(indir, 'rbindfill.r')
		infile1  = path.join(indir, 'rbindfill1_1.txt')
		infile2  = path.join(indir, 'rbindfill1_2.txt')
		outfile  = path.join(testdir, 'rbindfill1.txt')
		exptfile = path.join(outdir, 'rbindfill1.txt')
		yield rfile, infile1, infile2, outfile, exptfile
	
	def testRbindfill(self, rfile, infile1, infile2, outfile, exptfile):
		with helpers.captured_output() as (out, err):
			runcmd(['Rscript', rfile, infile1, infile2, outfile])
		self.assertFileEqual(outfile, exptfile)
	
if __name__ == '__main__':
	unittest.main(verbosity = 2)