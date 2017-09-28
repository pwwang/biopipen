import unittest, addPath

from os import path
from subprocess import Popen
from tempfile import gettempdir
from bioprocs.utils import plot

tmpdir = gettempdir()
class TestUtils (unittest.TestCase):

	def testPlotHeatmapR(self):
		scriptfile = path.join(tmpdir, 'plotHeatmap.r')
		outfile1    = path.join(tmpdir, 'plotHeatmap1.png')
		outfile2    = path.join(tmpdir, 'plotHeatmap2.png')
		outfile3    = path.join(tmpdir, 'plotHeatmap3.png')
		outfile4    = path.join(tmpdir, 'plotHeatmap4.png')
		outfile5    = path.join(tmpdir, 'plotHeatmap5.png')
		outfile6    = path.join(tmpdir, 'plotHeatmap6.png')
		outfile7    = path.join(tmpdir, 'plotHeatmap7.png')
		outfile8    = path.join(tmpdir, 'plotHeatmap8.png')
		outfile9    = path.join(tmpdir, 'plotHeatmap9.png')
		outfile10   = path.join(tmpdir, 'plotHeatmap10.png')
		outfile11   = path.join(tmpdir, 'plotHeatmap11.png')
		with open(scriptfile, 'w') as f:
			f.write(plot.heatmap.r)
			f.write("""
			m1 = matrix(rnorm(100), ncol = 10)
			m2 = matrix(rnorm(1000), ncol = 10)
			m3 = matrix(rnorm(1000), ncol = 100)
			colnames(m3) = paste("ThisIsAVeryLongColumnName", 1:100)
			m4 = matrix(rnorm(1000), ncol = 50)
			m51 = matrix(rnorm(50, 1, .1), ncol = 5)
			m52 = matrix(rnorm(70, 2, .1), ncol = 7)
			m5  = cbind(m51, m52)
			m5  = m5[, sample(12)]
			m61 = matrix(rnorm(30, 2, .1), ncol = 5)
			m6  = rbind(m51, m61)
			m6  = m6[sample(16), ]
			plotHeatmap(m1, "%s")
			plotHeatmap(m2, "%s")
			# add border
			plotHeatmap(m3, "%s", ggs = list(geom_tile(aes(fill=values), color="black")))
			plotHeatmap(m4, "%s")
			plotHeatmap(m5, "%s")
			plotHeatmap(m6, "%s")
			plotHeatmap(m1, "%s", ggs = list(theme(axis.text.x = element_blank())))
			plotHeatmap(m2, "%s", ggs = list(theme(axis.text.y = element_blank())))
			plotHeatmap(m3, "%s", ggs = list(theme(axis.text.x = element_blank(), axis.text.y = element_blank())))
			plotHeatmap(m4, "%s", ggs = list(theme(legend.position = "none")))
			plotHeatmap(m5, "%s", dendro = FALSE)
			""" % (outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, outfile11))
		rc = Popen(['Rscript', scriptfile]).wait()
		self.assertEqual(rc, 0)

	def testPlotBoxplot(self):
		scriptfile = path.join(tmpdir, 'plotBoxplot.r')
		outfile1    = path.join(tmpdir, 'plotBoxplot1.png')
		outfile2    = path.join(tmpdir, 'plotBoxplot2.png')
		with open(scriptfile, 'w') as f:
			f.write(plot.boxplot.r)
			f.write("""
			m1 = matrix(rnorm(100), ncol = 10)
			m2 = matrix(rnorm(1000), ncol = 10)
			plotBoxplot(m1, "%s")
			# ylab
			plotBoxplot(m2, "%s", ggs=list(ylab("Values"), theme(axis.title.y = element_text(angle = 90))))
			""" % (outfile1, outfile2))
		rc = Popen(['Rscript', scriptfile]).wait()
		self.assertEqual(rc, 0)

	def testPlotHist(self):
		scriptfile = path.join(tmpdir, 'plotHist.r')
		outfile1    = path.join(tmpdir, 'plotHist1.png')
		outfile2    = path.join(tmpdir, 'plotHist2.png')
		with open(scriptfile, 'w') as f:
			f.write(plot.hist.r)
			f.write("""
			m1 = matrix(rnorm(100), ncol = 10)
			m2 = matrix(rnorm(1000), ncol = 10)
			plotHist(m1, "%s")
			# ylab
			plotHist(m2, "%s", ggs=list(ylab("Values"), theme(axis.title.y = element_text(angle = 90))))
			""" % (outfile1, outfile2))
		rc = Popen(['Rscript', scriptfile]).wait()
		self.assertEqual(rc, 0)

	def testPlotVolPlot(self):
		scriptfile = path.join(tmpdir, 'plotVolplot.r')
		outfile1    = path.join(tmpdir, 'plotVolplot.png')
		with open(scriptfile, 'w') as f:
			f.write(plot.volplot.r)
			f.write("""
			logFC = rnorm(100, 0, 4)
			FDR   = abs(rnorm(100, .01, .005))
			logFCCut = rep(2, 100)
			FDRCut = rep(.01, 100)
			m = data.frame(logFC, FDR, logFCCut, FDRCut)
			plotVolplot(m, "%s")
			""" % (outfile1))
		rc = Popen(['Rscript', scriptfile]).wait()
		self.assertEqual(rc, 0)

if __name__ == '__main__':
	unittest.main(failfast=True, verbosity = 2)