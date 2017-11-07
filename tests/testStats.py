import addPath, unittest
import random
from os import path
from pyppl import PyPPL
from bioprocs import params
from bioprocs.stats import pMetaPval, pMetaPval1
from bioprocs.common import pFiles2Dir, pMergeFiles

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)
class testStats (unittest.TestCase):
	def test1MetaPval(self):
		rnames  = ['ROW'+str(i+1) for i in range(20)]
		infiles = []
		for i in range(10):
			infile = path.join(params.tmpdir.value, 'testMetaPval%s.txt' % i)
			infiles.append(infile)
			with open(infile, 'w') as fout:
				fout.write("ROW\tPVAL\n")
				for rname in random.sample(rnames, random.randint(10, 20)):
					fout.write("%s\t%s\n" % (rname, random.random()))
		pFiles2Dir.input = [infiles]
		pMetaPval.depends = pFiles2Dir
		PyPPL().start(pFiles2Dir).run()

	def test2MetaPval1(self):
		pMergeFiles.input = [[path.join(params.tmpdir.value, 'testMetaPval%s.txt' % i) for i in range(10)]]
		pMergeFiles.args.header = True

		pMetaPval1.depends = pMergeFiles
		
		PyPPL().start(pMergeFiles).run()

if __name__ == '__main__':
	unittest.main(failfast = True)