import addPath, unittest, random
from pyppl import PyPPL
from os import path
from tempfile import gettempdir
from bioprocs.plot import pHeatmap

tmpdir = gettempdir()
class testPlot (unittest.TestCase):
	
	def testpHeatmap(self):
		# input file
		infile = path.join(tmpdir, 'testpHeatmap.txt')
		rows   = ['ROW' + str(i+1) for i in range(100)]
		cols   = ['COL' + str(i+1) for i in range(20)]
		with open(infile, 'w') as f:
			f.write('\t'.join(cols) + '\n')
			for row in rows:
				f.write(row)
				for col in cols:
					f.write('\t' + str(random.randint(0, 100)))
				f.write('\n')
		pHeatmap.input = [infile]
		PyPPL().start(pHeatmap).run()

if __name__ == '__main__':
	unittest.main()
		