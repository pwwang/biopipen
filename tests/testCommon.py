import os, sys, unittest, addPath

from pyppl import pyppl, doct, channel, proc
import inspect
print inspect.getfile(pyppl)
from bioprocs.common import pSort, pFiles2List, pFiles2Dir

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestCommon (unittest.TestCase):

	def test (self):
		pFiles2List.input = {pFiles2List.input: [["c", "d", "e"], ["3", "4", "5"], ['1','1','1','1']]}
		
		pSort.depends = pFiles2List
		pSort.args.params = '-k1r'
		
		pChangeName = proc()
		pChangeName.depends = pSort
		pChangeName.input   = "infile:file"
		pChangeName.output  = "outfile:file:sameName.txt"
		pChangeName.script  = 'cp "{{infile}}" {{outfile}}'
		
		pFiles2Dir.depends = pChangeName
		pFiles2Dir.input   = {pFiles2Dir.input: lambda ch: [ch.toList()]}
		
		pyppl().starts(pFiles2List).run()

		
if __name__ == '__main__':
	unittest.main(failfast=True)