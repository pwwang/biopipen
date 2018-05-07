import helpers, testly, json
from os import path
from pyppl.templates import Template
from bioprocs.utils import runcmd

class TestUtilsPlot(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestUtilsPlot')

	def dataProvider_testScatter(self):
		datafile = path.join(self.indir, 'plot.scatter.data')
		plotfile1 = path.join(self.testdir, 'plot.scatter1.png')
		plotfile2 = path.join(self.testdir, 'plot.scatter2.png')
		plotfile3 = path.join(self.testdir, 'plot.scatter3.png')
		outfile1 = path.join(self.outdir, 'plot.scatter1.png')
		outfile3 = path.join(self.outdir, 'plot.scatter3.png')
		yield testly.Data(datafile, plotfile1, x = 1, y = 2, params = {}, outfile = outfile1)
		yield testly.Data(datafile, plotfile2, x = 'wt', y = 'mpg', params = {}, outfile = outfile1)
		yield testly.Data(datafile, plotfile3, x = 'wt', y = 'mpg', params = {
			'add': 'reg.line',
			'add.params': dict(color = "blue", fill = "lightgray"),
			'conf.int': True, # Add confidence interval
			'cor.coef': True, # Add correlation coefficient. see ?stat_cor
			'cor.coeff.args': {'method': "pearson", 'label.x': 3, 'label.sep': ", "}
		}, outfile = outfile3)

	def testScatter(self, datafile, plotfile, x, y, params, outfile):
		rscript = path.join(self.indir, 'plot.scatter.r')
		runcmd(['Rscript', rscript, datafile, plotfile, x, y, json.dumps(params)])
		self.assertFileEqual(plotfile, outfile)


if __name__ == '__main__':
	testly.main()
