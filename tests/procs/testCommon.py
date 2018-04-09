import helpers, unittest
from os import path
from glob import glob
from pyppl import PyPPL
from bioprocs.common import pSort, pFiles2Dir, pStr2File, pAddHeader, pMergeFiles, pFile2Proc, pHead, pTail, pPrepend, pAppend, pSplitRows
from helpers import config

class TestCommon (helpers.TestCase):
	
	def dataProvider_testSort(self, testdir, indir, outdir):
		infile = path.join(indir, 'sort.txt')
		yield 't1', infile, {}, {'k': 2}, True, False, path.join(outdir, 'sort.txt')
		yield 't2', infile, {'skip': 3}, {'k': 2}, True, False, path.join(outdir, 'sort-skip.txt')
		yield 't3', infile, {}, {'k': 2}, False, False, path.join(outdir, 'sort-case.txt')
		yield 't4', infile, {}, {'k': 1}, True, True, path.join(outdir, 'sort-unique.txt')
		yield 't5', infile, {'delimit': '|'}, {'k': 2}, True, True, path.join(outdir, 'sort-delimit.txt')
		
	def testSort(self, tag, infile, inopts, params, case, unique, exptfile):
		pSortTest = pSort.copy(tag = tag)
		pSortTest.input = [infile]
		pSortTest.args.inopts.update(inopts)
		pSortTest.args.params.update(params)
		pSortTest.args.case = case
		pSortTest.args.unique = unique
		PyPPL(config).start(pSortTest).run()
		self.assertFileEqual(pSortTest.channel.get(), exptfile)
		
	def dataProvider_testFiles2Dir(self, testdir, indir, outdir):
		files1 = glob(path.join(indir, 'files2dir*.txt'))
		yield 't1', files1, [path.basename(f) for f in files1]
		files2 = glob(path.join(indir, 'files2dirSameBn*', 'samefile.txt'))
		yield 't2', files2, ['samefile.txt' if i == 0 else 'samefile[%s].txt' % (i) for i, f in enumerate(files2)]
	
	def testFiles2Dir(self, tag, files, outfiles):
		pFiles2DirTest = pFiles2Dir.copy(tag = tag)
		pFiles2DirTest.input = [files]
		PyPPL(config).start(pFiles2DirTest).run()
		for outfile in outfiles:
			self.assertTrue(path.isfile(path.join(pFiles2DirTest.channel.get(), outfile)))
	
	def dataProvider_testStr2File(self, testdir, indir):
		infile = path.join(indir, 'sort.txt')
		yield 't1', infile
	
	def testStr2File(self, tag, infile):
		with open(infile) as f:
			instr = f.read()
		pStr2FileTest = pStr2File.copy(tag = tag)
		pStr2FileTest.input = [instr]
		PyPPL(config).start(pStr2FileTest).run()
		self.assertFileEqual(pStr2FileTest.channel.get(), infile)
	
	def dataProvider_testAddHeader(self, testdir, indir, outdir):
		infile1 = path.join(indir, 'sort.txt')
		infile2 = path.join(indir, 'sort.txt')
		outfile = path.join(outdir, 'addheader.txt')
		yield 't1', infile1, infile2, outfile
	
	def testAddHeader(self, tag, infile1, infile2, outfile):
		pAddHeaderTest = pAddHeader.copy(tag = tag)
		pAddHeaderTest.input = infile1, infile2
		PyPPL(config).start(pAddHeaderTest).run()
		self.assertFileEqual(pAddHeaderTest.channel.get(), outfile)
	
	def dataProvider_testMergeFiles(self, testdir, indir, outdir):
		files   = [path.join(indir, 'sort.txt')] * 2
		inopts  = {'skip': [0, 2], 'comment': ''}
		outfile = path.join(outdir, 'mergefiles.txt')
		yield 't1', files, inopts, 1000, outfile
		yield 't2', files, inopts, 1, outfile
	
	def testMergeFiles(self, tag, files, inopts, maxopen, outfile):
		pMergeFilesTest = pMergeFiles.copy(tag = tag)
		pMergeFilesTest.input = [files]
		pMergeFilesTest.args.maxopen = maxopen
		pMergeFilesTest.args.inopts.update(inopts)
		PyPPL(config).start(pMergeFilesTest).run()
		self.assertFileEqual(pMergeFilesTest.channel.get(), outfile)
		
	def dataProvider_testFile2Proc(self, testdir, indir, outdir):
		infile = path.join(indir, 'sort.txt')
		yield 't1', infile
	
	def testFile2Proc(self, tag, infile):
		pFile2ProcTest = pFile2Proc.copy(tag = tag)
		pFile2ProcTest.input = [infile]
		PyPPL(config).start(pFile2ProcTest).run()
		self.assertTrue(path.samefile(pFile2ProcTest.channel.get(), infile))
		
	def dataProvider_testHead(self, testdir, indir, outdir):
		infile1 = path.join(indir, 'sort.txt')
		outfile1 = path.join(outdir, 'head.txt')
		yield 't1', infile1, outfile1
	
	def testHead(self, tag, infile, outfile):
		pHeadTest = pHead.copy(tag = tag)
		pHeadTest.input = [infile]
		pHeadTest.args.n = 5
		PyPPL(config).start(pHeadTest).run()
		self.assertFileEqual(pHeadTest.channel.get(), outfile)
		
	def dataProvider_testTail(self, testdir, indir, outdir):
		infile1 = path.join(indir, 'sort.txt')
		outfile1 = path.join(outdir, 'tail.txt')
		yield 't1', infile1, outfile1
	
	def testTail(self, tag, infile, outfile):
		pTailTest = pTail.copy(tag = tag)
		pTailTest.input = [infile]
		pTailTest.args.n = '+5'
		PyPPL(config).start(pTailTest).run()
		self.assertFileEqual(pTailTest.channel.get(), outfile)
	
	def dataProvider_testPrepend(self, testdir, indir, outdir):
		infile = path.join(indir, 'sort.txt')
		outfile = path.join(outdir, 'prepend.txt')
		yield 't1', '#whatever to be prepended\n', infile, outfile
	
	def testPrepend(self, tag, instr, infile, outfile):
		pPrependTest = pPrepend.copy(tag = tag)
		pPrependTest.input = instr, infile
		PyPPL(config).start(pPrependTest).run()
		self.assertFileEqual(pPrependTest.channel.get(), outfile)
		
	def dataProvider_testAppend(self, testdir, indir, outdir):
		infile = path.join(indir, 'sort.txt')
		outfile = path.join(outdir, 'append.txt')
		yield 't1', '#whatever to be appended\n', infile, outfile
	
	def testAppend(self, tag, instr, infile, outfile):
		pAppendTest = pAppend.copy(tag = tag)
		pAppendTest.input = instr, infile
		PyPPL(config).start(pAppendTest).run()
		self.assertFileEqual(pAppendTest.channel.get(), outfile)
		
	def dataProvider_testSplitRows(self, testdir, indir, outdir):
		infile  = path.join(indir, 'sort.txt')
		outfile = path.join(outdir, 'splitrows.rows')
		args    = {'skip':0, 'cnames':True, 'n':8}
		yield 't1', infile, args, outfile
	
	def testSplitRows(self, tag, infile, args, outfile):
		pSplitRowsTest = pSplitRows.copy(tag = tag)
		pSplitRowsTest.input = [infile]
		pSplitRowsTest.args.update(args)
		PyPPL(config).start(pSplitRowsTest).run()
		self.assertDirEqual(pSplitRowsTest.channel.get(), outfile)
	
if __name__ == '__main__':
	unittest.main(failfast = True)