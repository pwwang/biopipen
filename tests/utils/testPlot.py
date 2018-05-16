import helpers, testly, json
from os import path
from pyppl.templates import Template
from bioprocs.utils import runcmd

class TestUtilsPlot(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestUtilsPlot')

	def dataProvider_testPoints(self):
		datafile = path.join(self.indir, 'plot.scatter.data')
		plotfile1 = path.join(self.testdir, 'plot.points1.png')
		plotfile2 = path.join(self.testdir, 'plot.points2.png')
		plotfile3 = path.join(self.testdir, 'plot.points3.png')
		outfile1 = path.join(self.outdir, 'plot.points1.png')
		outfile2 = path.join(self.outdir, 'plot.points2.png')
		outfile3 = path.join(self.outdir, 'plot.points3.png')
		yield testly.Data(datafile, plotfile1, x = 1, y = 2, params = {}, ggs = {}, outfile = outfile1)
		yield testly.Data(datafile, plotfile2, x = 1, y = 2, params = {}, ggs = {
			'geom_smooth': {}
		}, outfile = outfile2)
		yield testly.Data(datafile, plotfile3, x = 2, y = 1, params = {}, ggs = {
			'geom_smooth': {
				'method': 'lm', 'se': False
			},
			'annotate': {
				'geom': 'label',
				'label': '123456',
				'x': 20,
				'y': 1,
				'vjust': 0,
				'hjust': 0
			}
		}, outfile = outfile3)

	def testPoints(self, datafile, plotfile, x, y, params, ggs, outfile):
		rscript = path.join(self.indir, 'plot.points.r')
		runcmd(['Rscript', rscript, datafile, plotfile, x, y, json.dumps(params), json.dumps(ggs)])
		self.assertFileEqual(plotfile, outfile)

	def dataProvider_testBoxplot(self):
		datafile = path.join(self.indir, 'plot.scatter.data')
		plotfile1 = path.join(self.testdir, 'plot.boxplot1.png')
		outfile1 = path.join(self.outdir, 'plot.boxplot1.png')
		yield testly.Data(datafile, plotfile1, x=3, outfile = outfile1)

	def testBoxplot(self, datafile, plotfile, x = 1, y = 2, params = {}, ggs = {}, outfile = ''):
		rscript = path.join(self.indir, 'plot.boxplot.r')
		runcmd(['Rscript', rscript, datafile, plotfile, x, y, json.dumps(params), json.dumps(ggs)])
		self.assertFileEqual(plotfile, outfile)

	def dataProvider_testHeatmap(self):
		datafile = path.join(self.indir, 'plot.heatmap.data')
		plotfile1 = path.join(self.testdir, 'plot.heatmap1.png')
		outfile1 = path.join(self.outdir, 'plot.heatmap1.png')
		plotfile2 = path.join(self.testdir, 'plot.heatmap2.png')
		outfile2 = path.join(self.outdir, 'plot.heatmap2.png')
		plotfile3 = path.join(self.testdir, 'plot.heatmap3.png')
		outfile3 = path.join(self.outdir, 'plot.heatmap3.png')
		plotfile4 = path.join(self.testdir, 'plot.heatmap4.png')
		outfile4 = path.join(self.outdir, 'plot.heatmap4.png')
		yield testly.Data(datafile, plotfile1, outfile = outfile1)
		yield testly.Data(datafile, plotfile2, params = {'dendro': 'row'}, outfile = outfile2)
		yield testly.Data(datafile, plotfile3, params = {'dendro': 'col'}, outfile = outfile3)
		yield testly.Data(datafile, plotfile4, params = {'dendro': False}, outfile = outfile4)

	def testHeatmap(self, datafile, plotfile, params = {}, ggs = {}, outfile = ''):
		rscript = path.join(self.indir, 'plot.heatmap.r')
		runcmd(['Rscript', rscript, datafile, plotfile, json.dumps(params), json.dumps(ggs)])
		self.assertFileEqual(plotfile, outfile)

	def dataProvider_testFreqpoly(self):
		datafile = path.join(self.indir, 'plot.scatter.data')
		plotfile1 = path.join(self.testdir, 'plot.freqpoly1.png')
		outfile1 = path.join(self.outdir, 'plot.freqpoly1.png')
		yield testly.Data(datafile, plotfile1, x=3, outfile = outfile1)

	def testFreqpoly(self, datafile, plotfile, x = 1, params = {}, ggs = {}, outfile = ''):
		rscript = path.join(self.indir, 'plot.freqpoly.r')
		runcmd(['Rscript', rscript, datafile, plotfile, x, json.dumps(params), json.dumps(ggs)])
		self.assertFileEqual(plotfile, outfile)

	def dataProvider_testHisto(self):
		datafile = path.join(self.indir, 'plot.scatter.data')
		plotfile1 = path.join(self.testdir, 'plot.histo1.png')
		outfile1 = path.join(self.outdir, 'plot.histo1.png')
		yield testly.Data(datafile, plotfile1, x=3, outfile = outfile1)

	def testHisto(self, datafile, plotfile, x = 1, params = {}, ggs = {}, outfile = ''):
		rscript = path.join(self.indir, 'plot.histo.r')
		runcmd(['Rscript', rscript, datafile, plotfile, x, json.dumps(params), json.dumps(ggs)])
		self.assertFileEqual(plotfile, outfile)

if __name__ == '__main__':
	testly.main()
