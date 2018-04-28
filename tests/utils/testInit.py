import helpers, unittest
from collections import OrderedDict
from os import path
from helpers import testdirs
from subprocess import check_output
from bioprocs.utils import runcmd, cmdargs, RuncmdException, mem2, Mem2Exception


class TestInit(helpers.TestCase):

	testdir, indir, outdir = testdirs('TestInit')
	
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
	
	def dataProvider_testRuncmdR(self):
		srcfile = path.join(self.indir, 'runcmd1.r')
		yield srcfile, True
		srcfile = path.join(self.indir, 'runcmd2.r')
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
		
	def dataProvider_testCmdargsR(self):
		yield path.join(self.indir, 'cmdargs1.r'), ''
		yield path.join(self.indir, 'cmdargs2.r'), "-a 1.0 -c --de='de fg' -b 2.0"
		yield path.join(self.indir, 'cmdargs3.r'), "--a 1.0 --c --de='de fg' --b 2.0"
		yield path.join(self.indir, 'cmdargs4.r'), "--a 1.0 --c --de 'de fg' --b 2.0"
		
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
	
	def dataProvider_testMem2R(self):
		rfile = path.join(self.indir, 'mem2.r')
		yield rfile, '123k', 'auto', '123K'
		yield rfile, '100M', 'auto', '100M'
		yield rfile, '100M', 'K', '102400K'
		yield rfile, '100M', 'Java', '-Xms12800K -Xmx100M'
		yield rfile, '8G', 'Java', '-Xms1G -Xmx8G'
	
	def testMem2R(self, rfile, mem, unit, ret):
		out = check_output(['Rscript', rfile, mem, unit])
		self.assertEqual(out, ret)
		
	def dataProvider_testCbindfill(self):
		rfile    = path.join(self.indir, 'cbindfill.r')
		infile1  = path.join(self.indir, 'cbindfill1_1.txt')
		infile2  = path.join(self.indir, 'cbindfill1_2.txt')
		outfile  = path.join(self.testdir, 'cbindfill1.txt')
		exptfile = path.join(self.outdir, 'cbindfill1.txt')
		yield rfile, infile1, infile2, outfile, exptfile
		
	def testCbindfill(self, rfile, infile1, infile2, outfile, exptfile):
		with helpers.captured_output():
			runcmd(['Rscript', rfile, infile1, infile2, outfile])
		self.assertFileEqual(outfile, exptfile)
	
	def dataProvider_testRbindfill(self):
		rfile    = path.join(self.indir, 'rbindfill.r')
		infile1  = path.join(self.indir, 'rbindfill1_1.txt')
		infile2  = path.join(self.indir, 'rbindfill1_2.txt')
		outfile  = path.join(self.testdir, 'rbindfill1.txt')
		exptfile = path.join(self.outdir, 'rbindfill1.txt')
		yield rfile, infile1, infile2, outfile, exptfile
	
	def testRbindfill(self, rfile, infile1, infile2, outfile, exptfile):
		with helpers.captured_output() as (out, err):
			runcmd(['Rscript', rfile, infile1, infile2, outfile])
		self.assertFileEqual(outfile, exptfile)
	
if __name__ == '__main__':
	unittest.main(verbosity = 2)