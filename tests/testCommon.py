import os, sys, unittest, addPath

from tempfile import gettempdir
from pyppl import PyPPL, Channel, Proc
from bioprocs.common import pSort, pFiles2List, pFiles2Dir

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

tmpdir = gettempdir()
class TestCommon (unittest.TestCase):

	def test (self):
		pFiles2ListIn = [[tmpdir + "/c", tmpdir + "/d", tmpdir + "/e"], [tmpdir + "/3", tmpdir + "/4", tmpdir + "/5"], [tmpdir + "/1",tmpdir + "/1",tmpdir + "/1",tmpdir + "/1"]]
		for l in pFiles2ListIn:
			for i in l: open(i, 'w').close()

		pFiles2List.input = pFiles2ListIn
		
		pSort.depends = pFiles2List
		pSort.args.params = '-k1r'
		
		pChangeName = Proc()
		pChangeName.depends = pSort
		pChangeName.input   = "infile:file"
		pChangeName.output  = "outfile:file:sameName.txt"
		pChangeName.script  = 'cp "{{in.infile}}" {{out.outfile}}'
		
		pFiles2Dir.depends = pChangeName
		pFiles2Dir.input   = lambda ch: [ch.flatten()]
		
		PyPPL().start(pFiles2List).run()

		
if __name__ == '__main__':
	unittest.main(failfast=True)