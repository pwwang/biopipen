import os, sys, unittest, addPath
from os import path
from tempfile import gettempdir
from pyppl import PyPPL, Channel, Proc
from bioprocs.common import pSort, pFiles2List, pFiles2Dir, pStr2File

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
		
if __name__ == '__main__':
	unittest.main(failfast=True)