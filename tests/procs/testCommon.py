import unittest
from os import path
from pyppl import PyPPL
from bioprocs.common import pSort, pFiles2Dir, pStr2File, pAddHeader, pMergeFiles
from helpers import getfile, procOK, config

class testCommon (unittest.TestCase):
	
	def testSort(self):
		name                 = 'sort.txt'
		pSort1               = pSort.copy()
		pSort1.input         = [getfile(name)]
		pSort1.args.params.k = 2
		PyPPL(config).start(pSort1).run()
		procOK(pSort1, 'sort.txt', self)
	
	def testSortSkip(self):
		name                 = 'sort.txt'
		pSort2               = pSort.copy()
		pSort2.input         = [getfile(name)]
		pSort2.args.skip     = 3
		pSort2.args.params.k = 2
		PyPPL(config).start(pSort2).run()
		procOK(pSort2, 'sort-skip.txt', self)

	def testSortCase(self):
		name                 = 'sort.txt'
		pSort3               = pSort.copy()
		pSort3.input         = [getfile(name)]
		pSort3.args.case     = False
		pSort3.args.params.k = 2
		PyPPL(config).start(pSort3).run()
		procOK(pSort3, 'sort-case.txt', self)

	def testSortUnique(self):
		name                      = 'sort.txt'
		pSortUnique               = pSort.copy()
		pSortUnique.input         = [getfile(name)]
		pSortUnique.args.params.k = 1
		pSortUnique.args.unique   = True
		PyPPL(config).start(pSortUnique).run()
		procOK(pSortUnique, 'sort-unique.txt', self)
	
	def testSortDelimit(self):
		name                       = 'sort.txt'
		pSortDelimit               = pSort.copy()
		pSortDelimit.input         = [getfile(name)]
		pSortDelimit.args.params.k = 2
		pSortDelimit.args.delimit  = '|'
		PyPPL(config).start(pSortDelimit).run()
		procOK(pSortDelimit, 'sort-delimit.txt', self)

	def testSortNoeline(self):
		name                       = 'sort.txt'
		pSortNoeline               = pSort.copy()
		pSortNoeline.input         = [getfile(name)]
		pSortNoeline.args.params.k = 1
		pSortNoeline.args.noeline  = True
		PyPPL(config).start(pSortNoeline).run()
		procOK(pSortNoeline, 'sort-noeline.txt', self)

	def testFiles2Dir(self):
		names = ['files2dir' + str(i+1) + '.txt' for i in range(5)]
		pFiles2Dir1 = pFiles2Dir.copy()
		pFiles2Dir1.input = [[getfile(name) for name in names]]
		PyPPL(config).start(pFiles2Dir1).run()
		procOK(pFiles2Dir1, names[0][:-4] + '.dir', self)

	def testFile2DirSameBn(self):
		names                    = ['files2dirSameBn'+ str(i+1) +'/samefile.txt' for i in range(3)]
		pFiles2Dir2              = pFiles2Dir.copy()
		pFiles2Dir2.input        = [[getfile(name) for name in names]]
		PyPPL(config).start(pFiles2Dir2).run()
		procOK(pFiles2Dir2, 'samefile.dir', self)

	def testStr2File(self):
		with open(getfile('sort.txt')) as f:
			instr = f.read()
		pStr2File.input = [instr]
		PyPPL(config).start(pStr2File).run()
		procOK(pStr2File, 'str2file.txt', self)

	def testAddHeader(self):
		pAddHeader.input = (getfile('sort.txt'), getfile('sort.txt'))
		PyPPL(config).start(pAddHeader).run()
		procOK(pAddHeader, 'addheader.txt', self)

	def testMergeFiles(self):
		pMergeFiles.input = [[getfile('sort.txt'), getfile('sort.txt')]]
		pMergeFiles.args.skip = [0, 2]
		pMergeFiles.args.comment = ''
		PyPPL(config).start(pMergeFiles).run()
		procOK(pMergeFiles, 'mergefiles.txt', self)
	
if __name__ == '__main__':
	unittest.main(failfast = True)