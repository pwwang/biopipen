import unittest, addPath

from os import path
from subprocess import Popen
from tempfile import gettempdir
from bioprocs.utils import plot, txt

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

	def testTxtGroupPython(self):
		scriptfile = path.join(tmpdir, 'txtGroup.py')
		groupfileE = path.join(tmpdir, 'groupE.txt')
		groupfile2 = path.join(tmpdir, 'group2.txt')
		groupfile3 = path.join(tmpdir, 'group3.txt')
		with open(groupfileE, 'w') as f:
			f.write("a\n")
			f.write("a\tb\n")
		with open(groupfile2, 'w') as f:
			f.write("Sample\tGroup\n")
			f.write("sample1\tgroup1\n")
			f.write("sample2\tgroup1\n")
			f.write("sample3\tgroup1\n")
			f.write("sample4\tgroup2\n")
			f.write("sample5\tgroup2\n")
			f.write("sample6\tgroup2\n")
		with open(groupfile3, 'w') as f:
			f.write("Patient\tGroup\tBatch\n")
			f.write("sample1\tp1\tgroup1\tbatch1\n")
			f.write("sample2\tp1\tgroup1\tbatch1\n")
			f.write("sample3\tp2\tgroup1\tbatch1\n")
			f.write("sample4\tp2\tgroup2\tbatch2\n")
			f.write("sample5\tp3\tgroup2\tbatch2\n")
			f.write("sample6\tp3\tgroup2\tbatch2\n")
		with open(scriptfile, 'w') as f:
			f.write(txt.sampleinfo.python)
			f.write("""
r2info, r2cols = txtSampleinfo("%s")
assert r2info == {
	'sample1': {'Group': 'group1'},
	'sample2': {'Group': 'group1'},
	'sample3': {'Group': 'group1'},
	'sample4': {'Group': 'group2'},
	'sample5': {'Group': 'group2'},
	'sample6': {'Group': 'group2'},
}
assert r2cols == {
	'Sample': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
	'Group': ['group1', 'group1', 'group1', 'group2', 'group2', 'group2']
}

r3info, r3cols = txtSampleinfo("%s")
assert r3info == {
	'sample1': {'Group': 'group1', 'Patient': 'p1', 'Batch': 'batch1'},
	'sample2': {'Group': 'group1', 'Patient': 'p1', 'Batch': 'batch1'},
	'sample3': {'Group': 'group1', 'Patient': 'p2', 'Batch': 'batch1'},
	'sample4': {'Group': 'group2', 'Patient': 'p2', 'Batch': 'batch2'},
	'sample5': {'Group': 'group2', 'Patient': 'p3', 'Batch': 'batch2'},
	'sample6': {'Group': 'group2', 'Patient': 'p3', 'Batch': 'batch2'},
}
assert r3cols == {
	'Sample': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
	'Patient': ['p1', 'p1', 'p2', 'p2', 'p3', 'p3'],
	'Group': ['group1', 'group1', 'group1', 'group2', 'group2', 'group2'],
	'Batch': ['batch1', 'batch1', 'batch1', 'batch2', 'batch2', 'batch2'],
}

try:
	txtSampleinfo("%s")
	raise RuntimeError()
except ValueError:
	pass

""" % (groupfile2, groupfile3, groupfileE))

		rc = Popen(['python', scriptfile]).wait()
		self.assertEqual(rc, 0)

	def testTxtGroupR(self):
		scriptfile = path.join(tmpdir, 'txtGroup.r')
		groupfileE = path.join(tmpdir, 'groupE.txt')
		groupfile2 = path.join(tmpdir, 'group2.txt')
		groupfile3 = path.join(tmpdir, 'group3.txt')
		with open(groupfileE, 'w') as f:
			f.write("a\n")
			f.write("a\tb\n")
		with open(groupfile2, 'w') as f:
			f.write("Sample\tGroup\n")
			f.write("sample1\tgroup1\n")
			f.write("sample2\tgroup1\n")
			f.write("sample3\tgroup1\n")
			f.write("sample4\tgroup2\n")
			f.write("sample5\tgroup2\n")
			f.write("sample6\tgroup2\n")
		with open(groupfile3, 'w') as f:
			f.write("Patient\tGroup\tBatch\n")
			f.write("sample1\tp1\tgroup1\tbatch1\n")
			f.write("sample2\tp1\tgroup1\tbatch1\n")
			f.write("sample3\tp2\tgroup1\tbatch1\n")
			f.write("sample4\tp2\tgroup2\tbatch2\n")
			f.write("sample5\tp3\tgroup2\tbatch2\n")
			f.write("sample6\tp3\tgroup2\tbatch2\n")
		with open(scriptfile, 'w') as f:
			f.write(txt.sampleinfo.r)
			f.write("""
				r2 = txtSampleinfo("%s")
				
				stopifnot( rownames(r2) == c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6') )
				stopifnot( r2$Group == c('group1', 'group1', 'group1', 'group2', 'group2', 'group2') )

				r3 = txtSampleinfo("%s")
				stopifnot( rownames(r3) == c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6') )
				stopifnot( r3$Group == c('group1', 'group1', 'group1', 'group2', 'group2', 'group2') )
				stopifnot( r3$Patient == c('p1', 'p1', 'p2', 'p2', 'p3', 'p3') )
				stopifnot( r3$Batch == c('batch1', 'batch1', 'batch1', 'batch2', 'batch2', 'batch2') )

				tryCatch({
					txtSampleinfo("%s")
					stop('Not happen')
				}, error = function(e){
					stopifnot(!'Not happen' %%in%% e)
				})
			""" % (groupfile2, groupfile3, groupfileE))

		rc = Popen(['Rscript', scriptfile]).wait()
		self.assertEqual(rc, 0)

if __name__ == '__main__':
	unittest.main(failfast=True, verbosity = 2)