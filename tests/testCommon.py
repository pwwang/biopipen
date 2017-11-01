import os, sys, unittest, addPath
from os import path
from tempfile import gettempdir
from pyppl import PyPPL, Channel, Proc
from bioprocs.common import pSort, pFiles2List, pFiles2Dir, pStr2File, pMergeFiles

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

tmpdir = gettempdir()
class TestCommon (unittest.TestCase):

	def test (self):
		pFiles2ListIn = [[tmpdir + "/c", tmpdir + "/d", tmpdir + "/e"], [tmpdir + "/3", tmpdir + "/4", tmpdir + "/5"], [tmpdir + "/1",tmpdir + "/1",tmpdir + "/1",tmpdir + "/1"]]
		for l in pFiles2ListIn:
			for i in l: 
				if not path.exists(i):
					open(i, 'w').close()

		pFiles2List.input = pFiles2ListIn
		pFiles2List.forks = 3

		pSort.depends         = pFiles2List
		pSort.forks           = 3
		pSort.args.skip       = 1
		pSort.args.params.k1r = True

		
		pChangeName = Proc()
		pChangeName.depends = pSort
		pChangeName.input   = "infile:file"
		pChangeName.output  = "outfile:file:sameName.txt"
		pChangeName.forks   = 3
		pChangeName.script  = 'cp "{{in.infile}}" {{out.outfile}}'
		
		pFiles2Dir.depends = pChangeName
		pFiles2Dir.input   = lambda ch: [ch.flatten()]
		
		PyPPL().start(pFiles2List).run()

	def testStr2File(self):
		pStr2File.input = ["a'b", "a, b"]
		pStr2File.forks = 2
		PyPPL().start(pStr2File).run()
		with open(pStr2File.channel.get()) as f, open(pStr2File.channel.get(1)) as f1:
			self.assertEqual(f.read().strip(), "a'b")
			self.assertEqual(f1.read().strip(), "a\nb")

	def testMergeFiles(self):
		file1              = path.join(tmpdir, 'mergefile1.txt')
		file2              = path.join(tmpdir, 'mergefile2.txt')
		file3              = path.join(tmpdir, 'mergefile3.txt')
		with open(file1, 'w') as f1, open(file2, 'w') as f2, open(file3, 'w') as f3:
			f1.write('# skiped\n')
			f1.write('# skiped2\n')
			f1.write('m1	1\n')
			f2.write('M	N\n')
			f2.write('m2	2\n')
			f3.write('m3	3\n')
		pMergeFiles.input       = [[file1, file2, file3]]
		pMergeFiles.args.header = [0,1,0]
		PyPPL().start(pMergeFiles).run()
		with open(pMergeFiles.channel.get()) as f:
			self.assertEqual(f.read().splitlines(), ['M	N', 'm1	1', 'm2	2', 'm3	3'])

		
if __name__ == '__main__':
	unittest.main(failfast=True)