import unittest
from os import path
from pyppl import PyPPL
from bioprocs.utils import plot, txt, helpers, mem2, runcmd, polling, genenorm
from helpers import getfile, procOK, config, utilTest

# utilTest(input, script, name, tplenvs, test, args = None)
class TestUtils (unittest.TestCase):

	def testPlotHeatmapR(self):
		utilTest(
			{'infile:file': getfile('heatmap1.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap1.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self
		)

		utilTest(
			{'infile:file': getfile('heatmap2.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap2.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self
		)

		utilTest(
			{'infile:file': getfile('heatmap3.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap3.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'ggs': ['r:geom_tile(aes(fill=values), color="black")']}
		)

		utilTest(
			{'infile:file': getfile('heatmap4.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap4.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self
		)

		utilTest(
			{'infile:file': getfile('heatmap5.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap5.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self
		)

		utilTest(
			{'infile:file': getfile('heatmap6.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap6.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self
		)

		utilTest(
			{'infile:file': getfile('heatmap1.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap7.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'ggs': ['r:theme(axis.text.x = element_blank())']}
		)

		utilTest(
			{'infile:file': getfile('heatmap2.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap8.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'ggs': ['r:theme(axis.text.y = element_blank())']}
		)

		utilTest(
			{'infile:file': getfile('heatmap3.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap9.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'ggs': ['r:theme(axis.text.x = element_blank(), axis.text.y = element_blank())']}
		)

		utilTest(
			{'infile:file': getfile('heatmap4.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap10.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'ggs': ['r:theme(legend.position = "none")']}
		)

		utilTest(
			{'infile:file': getfile('heatmap5.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap11.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'dendro': {'dendro': False}}
		)

		utilTest(
			{'infile:file': getfile('heatmap1.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap12.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'dendro': {'dendro': 'col', 'rows': ['2', '3', '4', '5', '6']}}
		)

		utilTest(
			{'infile:file': getfile('heatmap1.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap13.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'dendro': {'dendro': 'row', 'cols': ['V2', 'V3', 'V4', 'V5', 'V6', 'V7']}}
		)

		utilTest(
			{'infile:file': getfile('heatmap1.txt')},
			getfile('plotHeatmap1.r'), 
			'plotHeatmap14.png', 
			{'plotHeatmap': plot.heatmap.r}, 
			self,
			{'dendro': {'dendro': False, 'rows': ['6', '7', '8', '9', '10'], 'cols': ['V7', 'V8', 'V9', 'V10']}}
		)

	def testTxtFilter(self):
		utilTest(
			{'infile:file': getfile('txtfilter.txt')},
			getfile('txtfilter.py'), 
			'txtfilter.txt', 
			{'txtFilter': txt.filter.py}, 
			self
		)
	
	def testTxtTransform(self):
		utilTest(
			{'infile:file': getfile('txtfilter.txt')},
			getfile('txttransform.py'), 
			'txttransform.txt', 
			{'txtTransform': txt.transform.py}, 
			self
		)

	def testRuncmdR(self):
		utilTest(
			{'in': [1]},
			getfile('runcmd.r'),
			'runcmdR.txt',
			{'runcmd': runcmd.r},
			self
		)
	
	def testPolling1stR(self):
		utilTest(
			{'in': [1, 2]},
			getfile('pollinglast.r'),
			'pollinglastR.txt',
			{'pollingLast': polling.last.r},
			self
		)

if __name__ == '__main__':
	unittest.main()#failfast=True)