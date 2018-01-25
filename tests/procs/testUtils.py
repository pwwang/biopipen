import unittest
from os import path
from pyppl import PyPPL
from bioprocs.utils import plot, txt, helpers, mem2, runcmd, polling, genenorm, read, write
from helpers import getfile, procOK, config, utilTest, cmdOK, fileOK
from bioprocs import params

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
	
	def testReadMetaPy(self):
		cmdOK([getfile('testReadMetaPy.py')], self)
	
	def testReadRecordPy(self):
		cmdOK([getfile('testReadRecordPy.py')], self)

	def testReadBasePy(self):
		cmdOK([getfile('testReadBasePy.py'), getfile('testReadBasePy.txt')], self)

	def testWriteBasePy(self):
		cmdOK([
			getfile('testWriteBasePy.py'), 
			getfile('testReadBasePy.txt'), 
			path.join(params.tmpdir.value, 'testWriteBasePy.txt')
		], self)
		fileOK(
			getfile('testWriteBasePy.txt', False), 
			path.join(params.tmpdir.value, 'testWriteBasePy.txt'), self)

	def testReadBed12Py(self):
		cmdOK([getfile('testReadBed12Py.py'), getfile('testReadBed12Py.bed12')], self)

	def testReadBedpePy(self):
		cmdOK([getfile('testReadBedpePy.py'), getfile('testReadBedpePy.bedpe')], self)

	def testReadBedxPy(self):
		print ' '.join([getfile('testReadBedxPy.py'), getfile('testReadBedxPy.bedx')])
		cmdOK([getfile('testReadBedxPy.py'), getfile('testReadBedxPy.bedx')], self)

	def testReadBedPy(self):
		cmdOK([getfile('testReadBedPy.py'), getfile('testReadBedxPy.bedx')], self)


	def testWriteBedxPy(self):
		cmdOK([
			getfile('testWriteBedxPy.py'), 
			getfile('testReadBedxPy.bedx'), 
			path.join(params.tmpdir.value, 'testWriteBedxPy.txt')
		], self)
		fileOK(
			getfile('testWriteBedxPy.txt', False), 
			path.join(params.tmpdir.value, 'testWriteBedxPy.txt'), self)		


if __name__ == '__main__':
	unittest.main()#failfast=True)