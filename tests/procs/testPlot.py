import unittest, helpers, testly
from os import path
from glob import glob
from pyppl import PyPPL, Box
from helpers import testdirs, config
from bioprocs.plot import pScatter, pHeatmap, pBoxplot, pHisto, pFreqpoly, pScatterCompare, pROC, pVenn, pPie

class TestPlot (helpers.TestCase):

	testdir, indir, outdir = testdirs('TestPlot')

	def dataProvider_testScatter(self):
		infile = path.join(self.indir, 'scattercompare.txt')
		outfile = path.join(self.outdir, 'scatter1.png')
		outfile2 = path.join(self.outdir, 'scatter2.png')
		yield 't1', infile, outfile, 'Method1', 'Method2', True, True
		yield 't2', infile, outfile2, 1, 2, True, True, {}, {
			'geom_smooth': {'method': 'lm'},
			'geom_label': {'label': 'r:corr', 'x': .5, 'y': .2}
		}, 'corr = round(cor(data[,1], data[,2]), 3)'

	def testScatter(self, tag, infile, outfile, x, y, cnames, rnames, params = {}, ggs = {}, helper = '', devpars = Box(res = 48, height = 300, width =300)):
		pScatterTest              = pScatter.copy(tag = tag)
		pScatterTest.input        = [infile]
		pScatterTest.args.cnames  = cnames
		pScatterTest.args.rnames  = rnames
		pScatterTest.args.x       = x
		pScatterTest.args.y       = y
		pScatterTest.args.devpars = devpars
		pScatterTest.args.params  = params
		pScatterTest.args.ggs     = ggs
		pScatterTest.args.helper  = helper

		PyPPL(config).start(pScatterTest).run()
		self.assertFileEqual(pScatterTest.channel.get(), outfile)

	def dataProvider_testpHeatmap(self):
		datafile = path.join(self.indir, 'heatmap.txt')
		outfile  = path.join(self.outdir, 'heatmap.png')
		yield testly.Data(
			tag      = 't1',
			datafile = datafile,
			outfile  = outfile
		)

	def testpHeatmap(self, tag, datafile, params = {}, ggs = {}, outfile = '', devpars = Box(res = 48, height = 300, width =300)):
		pHeatmapTest              = pHeatmap.copy(tag = tag)
		pHeatmapTest.input        = [datafile]
		pHeatmapTest.args.params  = params
		pHeatmapTest.args.devpars = devpars
		pHeatmapTest.args.ggs     = ggs
		PyPPL(config).start(pHeatmapTest).run()
		self.assertFileEqual(pHeatmapTest.channel.get(), outfile)

	def dataProvider_testpBoxplot(self):
		datafile = path.join(self.indir, 'boxplot.txt')
		outfile  = path.join(self.outdir, 'boxplot.png')
		yield testly.Data(
			tag      = 't1',
			x        = 3,
			y        = 2,
			datafile = datafile,
			outfile  = outfile
		)

	def testpBoxplot(self, tag, datafile, x, y, params = {}, ggs = {}, outfile = '', devpars = Box(res = 48, height = 300, width =300)):
		pBoxplotTest              = pBoxplot.copy(tag = tag)
		pBoxplotTest.input        = [datafile]
		pBoxplotTest.args.params  = params
		pBoxplotTest.args.devpars = devpars
		pBoxplotTest.args.x       = x
		pBoxplotTest.args.y       = y
		pBoxplotTest.args.ggs     = ggs
		PyPPL(config).start(pBoxplotTest).run()
		self.assertFileEqual(pBoxplotTest.channel.get(), outfile)

	def dataProvider_testpHisto(self):
		datafile = path.join(self.indir, 'boxplot.txt')
		outfile  = path.join(self.outdir, 'histo.png')
		yield testly.Data(
			tag      = 't1',
			x        = 2,
			datafile = datafile,
			outfile  = outfile
		)

	def testpHisto(self, tag, datafile, x, params = {}, ggs = {}, outfile = '', devpars = Box(res = 48, height = 300, width =300)):
		pHistoTest              = pHisto.copy(tag = tag)
		pHistoTest.input        = [datafile]
		pHistoTest.args.params  = params
		pHistoTest.args.devpars = devpars
		pHistoTest.args.x       = x
		pHistoTest.args.ggs     = ggs
		PyPPL(config).start(pHistoTest).run()
		self.assertFileEqual(pHistoTest.channel.get(), outfile)

	def dataProvider_testpFreqpoly(self):
		datafile = path.join(self.indir, 'boxplot.txt')
		outfile  = path.join(self.outdir, 'freqpoly.png')
		yield testly.Data(
			tag      = 't1',
			x        = 3,
			datafile = datafile,
			outfile  = outfile
		)

	def testpFreqpoly(self, tag, datafile, x, params = {}, ggs = {}, outfile = '', devpars = Box(res = 48, height = 300, width =300)):
		pFreqpolyTest              = pFreqpoly.copy(tag = tag)
		pFreqpolyTest.input        = [datafile]
		pFreqpolyTest.args.params  = params
		pFreqpolyTest.args.devpars = devpars
		pFreqpolyTest.args.x       = x
		pFreqpolyTest.args.ggs     = ggs
		PyPPL(config).start(pFreqpolyTest).run()
		self.assertFileEqual(pFreqpolyTest.channel.get(), outfile)

	def dataProvider_testpScatterCompare(self):
		datafile = path.join(self.indir, 'scattercompare.txt')
		outfile  = path.join(self.outdir, 'scattercompare.png')
		yield testly.Data(
			tag = 't1',
			datafile = datafile,
			outfile = outfile
		)

	def testpScatterCompare(self, tag, datafile, x = 1, y = 2, params = {}, ggs = {}, outfile = '', devpars = Box(res = 48, height = 300, width =300)):
		pScatterCompareTest = pScatterCompare.copy(tag = tag)
		pScatterCompareTest.input = [datafile]
		pScatterCompareTest.args.x = x
		pScatterCompareTest.args.y = y
		pScatterCompareTest.args.params = params
		pScatterCompareTest.args.ggs = ggs
		pScatterCompareTest.args.devpars = devpars
		PyPPL(config).start(pScatterCompareTest).run()
		self.assertFileEqual(pScatterCompareTest.channel.get(), outfile)

	def dataProvider_testpROC(self):
		datafile1 = path.join(self.indir, 'roc-single.txt')
		outauc1   = path.join(self.outdir, 'roc-single/auc.txt')
		outpng1   = path.join(self.outdir, 'roc-single/roc.png')
		datafile2 = path.join(self.indir, 'roc-multi.txt')
		outauc2   = path.join(self.outdir, 'roc-multi-combine/auc.txt')
		outpng2   = path.join(self.outdir, 'roc-multi-combine/roc.png')
		datafile3 = datafile2
		outauc3   = path.join(self.outdir, 'roc-multi-nocombine/auc.txt')
		outpng3   = path.join(self.outdir, 'roc-multi-nocombine/Method1.roc.png')
		yield 'single', datafile1, outauc1, outpng1, {
			'cnames': False,
			'rnames': False
		}
		yield 'mcombine', datafile2, outauc2, outpng2, {
			'cnames': True,
			'rnames': False,
			'params': {'combine': True}
		}
		yield 'nocombine', datafile3, outauc3, outpng3, {
			'cnames': True,
			'rnames': False,
			'params': {'combine': False},
		}

	def testpROC(self, tag, datafile, outauc, outpng, args = None):
		pROCTest = pROC.copy(tag = tag)
		args = args or {}
		if 'params' not in args: args['params'] = {}
		params = {}
		params.update(pROCTest.args.params)
		params.update(args['params'])
		pROCTest.args.update(args)
		pROCTest.args.params.update(params)
		pROCTest.args.devpars = {'res': 96, 'width': 200, 'height': 200}
		pROCTest.input = [datafile]
		PyPPL(config).start(pROCTest).run()
		self.assertFileEqual(glob(path.join(pROCTest.channel.get(), '*.txt'))[0], outauc)
		self.assertFileEqual(glob(path.join(pROCTest.channel.get(), '*.png'))[0], outpng)

	def dataProvider_testpVenn(self):
		datafile1 = path.join(self.indir, 'venn-venn.txt')
		datafile2 = path.join(self.indir, 'venn-upset.txt')
		outfile1  = path.join(self.outdir, 'venn-venn.venn.png')
		outfile2  = path.join(self.outdir, 'venn-upset.venn.png')
		args1     = {'rnames': True}
		args2     = args1
		yield 'venn', datafile1, outfile1, args1
		yield 'upset', datafile2, outfile2, args2

	def testpVenn(self, tag, datafile, outfile, args):
		pVennTest = pVenn.copy(tag = tag)
		pVennTest.input = [datafile]
		pVennTest.args.update(args)
		pVennTest.args.devpars = {'res': 96, 'width': 200, 'height': 200}
		PyPPL(config).start(pVennTest).run()
		self.assertFileEqual(pVennTest.channel.get(), outfile)

	def dataProvider_testpPie(self):
		datafile1 = path.join(self.indir, 'pie-direct-num.txt')
		datafile2 = path.join(self.indir, 'pie-item-presence.txt')
		args1 = {}
		args2 = {'rnames': True}
		outfile1 = path.join(self.outdir, 'pie-direct-num.pie.png')
		outfile2 = path.join(self.outdir, 'pie-item-presence.pie.png')
		yield 't1', datafile1, args1, outfile1
		yield 't2', datafile2, args2, outfile2

	def testpPie(self, tag, datafile, args, outfile):
		pPieTest = pPie.copy(tag = tag)
		pPieTest.input = [datafile]
		pPieTest.args.update(args)
		pPieTest.args.devpars = {'res': 96, 'width': 200, 'height': 200}
		PyPPL(config).start(pPieTest).run()
		self.assertFileEqual(pPieTest.channel.get(), outfile)


if __name__ == '__main__':
	testly.main(failfast = True)
