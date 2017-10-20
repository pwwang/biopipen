import os, sys, unittest, addPath
from os import path
from pyppl import PyPPL, Channel, Proc
from bioprocs.seq import pPromoters
from bioprocs.common import pStr2File

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestSeq (unittest.TestCase):

	def testpPromoters(self):
		pStr2File.input = ['TP53']
		pPromoters.depends = pStr2File
		PyPPL().start(pStr2File).run()

		
if __name__ == '__main__':
	unittest.main(failfast=True)