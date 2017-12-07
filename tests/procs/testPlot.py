import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.plot import pHeatmap, pScatterCompare, pROC

class testPlot (unittest.TestCase):
	
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

	# TODO:
	# - boxplot
	# - scatter plot
	# - venn plot

if __name__ == '__main__':
	unittest.main(failfast = True)
		