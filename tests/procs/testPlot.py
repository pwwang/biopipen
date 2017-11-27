import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.plot import pHeatmap, pScatterCompare

class testPlot (unittest.TestCase):
	
	def testpHeatmap(self):
		pHeatmap.input = [getfile('heatmap.txt')]
		PyPPL(config).start(pHeatmap).run()
		procOK(pHeatmap, 'heatmap.png', self)

	def testpScatterCompare(self):
		pScatterCompare.input    = [getfile('scattercompare.txt')]
		PyPPL(config).start(pScatterCompare).run()
		procOK(pScatterCompare, 'scattercompare.png', self)

	# TODO:
	# - boxplot
	# - scatter plot
	# - venn plot

if __name__ == '__main__':
	unittest.main()
		