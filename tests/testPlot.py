import addPath, unittest, random
from pyppl import PyPPL
from os import path
from tempfile import gettempdir
from bioprocs.plot import pHeatmap, pScatterCompare

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

	def testpScatterCompare(self):
		# input file
		infile = path.join(tmpdir, 'testpScatterCompare.txt')
		rows   = ['ROW' + str(i+1) for i in range(100)]
		cols   = ['Method1', 'Method2', 'Group']
		with open(infile, 'w') as f:
			f.write('\t'.join(cols) + '\n')
			for row in rows:
				r2 = random.random()
				r1 = random.random() 
				r1 = r1 if r1 > r2 else random.random()
				r1 = r1 if r1 > r2 else random.random()
				f.write('%s\t%s\t%s\t%s\n' % (row, r1, r2, random.choice([0, 1])))
		pScatterCompare.input    = [infile]
		#pScatterCompare.args.ggs = ['r:geom_point(aes(color = factor(Group)))']
		PyPPL().start(pScatterCompare).run()
		


if __name__ == '__main__':
	unittest.main()
		