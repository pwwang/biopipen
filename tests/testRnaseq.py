import unittest, addPath

import random
from os import makedirs, path
from tempfile import gettempdir
from pyppl import PyPPL, Channel, Box, Proc
from bioprocs.rnaseq import pExpdir2Matrix, pBatchEffect, pRawCounts2, pDeg

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)


class TestRnaseq (unittest.TestCase):

	wdir = path.join(path.dirname(path.dirname(__file__)), 'workdir')
	data = {}

	@unittest.skipIf(path.exists(path.join(path.dirname(path.dirname(__file__)), 'workdir', 'data')), 'Data already generated.')
	def test0genData(self):
		nb1     = 10
		nb2     = 10
		genes   = ['Gene' + str(i+1) for i in range(1000)]
		genes.extend(['Sample', 'Composite', '__unmapped__'])
		datadir = path.join(self.wdir, 'data')
		batfile = path.join(self.wdir, 'data', 'batch.txt')
		if not path.exists(datadir): makedirs(datadir)
		bf = open(batfile, 'w') 
		for i, si in enumerate([range(nb1), range(nb2)]):
			for s in si:
				random.seed(100 * i + s)
				datafile = path.join(datadir, 'Sample' + str(i*10 + s + 1) + '.expr')
				bf.write('%s\t%d\n' % ('Sample' + str(i*10 + s + 1), i + 1))
				with open(datafile, 'w') as df:
					for gene in genes:
						df.write('%s\t%.3f\n' % (gene, random.normalvariate(10 - i * 3, 7)))
		bf.close()

	def test1pExpdir2Matrix(self):
		pExpdir2Matrix.input = [path.join(self.wdir, 'data')]
		pExpdir2Matrix.expect = "grep Gene100 {{out.outfile}} && grep Sample20 {{out.outfile}} && !(grep Composite {{out.outfile}})"
		pExpdir2Matrix.args.pattern  = '*.expr'
		pExpdir2Matrix.args.boxplot  = True
		pExpdir2Matrix.args.heatmap  = True
		pExpdir2Matrix.args.histplot = True
		PyPPL().start(pExpdir2Matrix).run()
		self.data['expfile'] = pExpdir2Matrix.channel.outfile.flatten()[0]

	def test2pBatchEffect(self):
		pBatchEffect.input = [(self.data['expfile'], path.join(self.wdir, 'data', 'batch.txt'))]
		pBatchEffect.args.boxplot  = True
		pBatchEffect.args.heatmap  = True
		pBatchEffect.args.histplot = True
		PyPPL().start(pBatchEffect).run()

	def test3pRawCounts2(self):
		pRawCounts2.input = [self.data['expfile']]
		pRawCounts2.args.boxplot  = True
		pRawCounts2.args.heatmap  = True
		pRawCounts2.args.histplot = True
		PyPPL().start(pRawCounts2).run()
		self.data['expfile'] = pExpdir2Matrix.channel.outfile.flatten()[0]

	def test4pDeg(self):
		pPos = Proc(desc = 'Make values positive.')
		pPos.input  = {"infile:file": self.data['expfile']}
		pPos.output = "outfile:file:{{in.infile|fn}}"
		pPos.script = "sed 's/-//g' {{in.infile}} > {{out.outfile}}"

		pDegPaired        = pDeg.copy()
		pDeg.depends      = pPos
		pDeg.args.heatmap = True
		pDeg.args.maplot  = True
		pDeg.input        = lambda ch: ch.cbind([path.join(self.wdir, 'data', 'batch.txt')])

		pDegPaired.depends      = pPos
		pDegPaired.args.heatmap = True
		pDegPaired.args.paired  = True
		pDegPaired.input        = lambda ch: ch.cbind([path.join(self.wdir, 'data', 'batch.txt')])
		PyPPL().start(pPos).run()


if __name__ == '__main__':
	unittest.main(failfast=True)