import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.plot import pHeatmap, pScatterCompare, pROC, pVenn, pPie

class TestPlot (unittest.TestCase):
	
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
		

	# TODO:
	# - boxplot
	# - scatter plot

if __name__ == '__main__':
	unittest.main(failfast = True)
		