import unittest, helpers, testly
from os import path
from pyppl import PyPPL, Box
from helpers import testdirs, config
from bioprocs.plot import pScatter
#from bioprocs.plot import pHeatmap, pScatterCompare, pROC, pVenn, pPie

class TestPlot (helpers.TestCase):

	testdir, indir, outdir = testdirs('TestPlot')

	def dataProvider_testScatter(self):
		infile = path.join(self.indir, 'scattercompare.txt')
		outfile = path.join(self.outdir, 'scatter1.png')
		outfile2 = path.join(self.outdir, 'scatter2.png')
		yield 't1', infile, outfile, 'Method1', 'Method2', True, True
		yield 't2', infile, outfile2, 1, 2, True, True, {
			'add': 'reg.line',
			'add.params': dict(color = "blue", fill = "lightgray"),
			'conf.int': True, # Add confidence interval
			'cor.coef': True, # Add correlation coefficient. see ?stat_cor
			'cor.coeff.args': {'method': "pearson", 'label.sep': ", "}
		}

	def testScatter(self, tag, infile, outfile, x, y, cnames, rnames, params = {}, devpars = Box(res = 48, height = 300, width =300)):
		pScatterTest = pScatter.copy(tag = tag)
		pScatterTest.input = [infile]
		pScatterTest.args.cnames = cnames
		pScatterTest.args.rnames = rnames
		pScatterTest.args.x = x
		pScatterTest.args.y = y
		pScatterTest.args.devpars = devpars
		pScatterTest.args.params = params
		PyPPL(config).start(pScatterTest).run()
		self.assertFileEqual(pScatterTest.channel.get(), outfile)

	'''
	def testpHeatmap(self):
		pHeatmap.input = [getfile('heatmap.txt')]
		PyPPL(config).start(pHeatmap).run()
		procOK(pHeatmap, 'heatmap.png', self)

	def testpScatterCompare(self):
		pScatterCompare.input    = [getfile('scattercompare.txt')]
		PyPPL(config).start(pScatterCompare).run()
		procOK(pScatterCompare, 'scattercompare.png', self)

	def testpROCsingle(self):
		pROC1 = pROC.copy()
		pROC1.input = [getfile('roc-single.txt')]
		pROC1.args.cnames = False
		pROC1.args.rnames = False
		PyPPL(config).start(pROC1).run()
		procOK(pROC1, 'roc-single', self)

	def testpROCmulti_combine(self):
		pROC2              = pROC.copy()
		pROC2.input        = [getfile('roc-multi.txt')]
		pROC2.args.cnames  = True
		pROC2.args.rnames  = False
		pROC2.args.combine = True
		PyPPL(config).start(pROC2).run()
		procOK(pROC2, 'roc-multi-combine', self)

	def testpROCmulti_nocombine(self):
		pROC3              = pROC.copy()
		pROC3.input        = [getfile('roc-multi.txt')]
		pROC3.args.cnames  = True
		pROC3.args.rnames  = False
		pROC3.args.combine = False
		PyPPL(config).start(pROC3).run()
		procOK(pROC3, 'roc-multi-nocombine', self)

	def testVennVenn(self):
		pVennVenn             = pVenn.copy()
		pVennVenn.input       = [getfile('venn-venn.txt')]
		pVennVenn.args.rnames = True
		PyPPL(config).start(pVennVenn).run()
		procOK(pVennVenn, 'venn-venn.venn.png', self)

	def testVennUpset(self):
		pVennUpset             = pVenn.copy()
		pVennUpset.input       = [getfile('venn-upset.txt')]
		pVennUpset.args.rnames = True
		PyPPL(config).start(pVennUpset).run()
		procOK(pVennUpset, 'venn-upset.venn.png', self)

	def testPieDirectNumbers(self):
		pPie1 = pPie.copy()
		pPie1.input = [getfile('pie-direct-num.txt')]
		PyPPL(config).start(pPie1).run()
		procOK(pPie1, 'pie-direct-num.pie.png', self)

	def testPieItemPresence(self):
		pPie2 = pPie.copy()
		pPie2.input = [getfile('pie-item-presence.txt')]
		pPie2.args.rnames = True
		PyPPL(config).start(pPie2).run()
		procOK(pPie2, 'pie-item-presence.pie.png', self)
	'''

	# TODO:
	# - boxplot
	# - scatter plot

if __name__ == '__main__':
	testly.main(failfast = True)
