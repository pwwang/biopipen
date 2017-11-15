import unittest, addPath

import random
from os import path, devnull, remove, makedirs
from subprocess import Popen, check_output
from tempfile import gettempdir
from bioprocs.utils import plot, txt, helpers, mem2, runcmd, polling, genenorm

tmpdir = gettempdir()
class TestUtils (unittest.TestCase):

	def setUp(self):
		self.devnull = open(devnull, 'w')

	def tearDown(self):
		if self.devnull:
			self.devnull.close()

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
		outfile12   = path.join(tmpdir, 'plotHeatmap12.png')
		outfile13   = path.join(tmpdir, 'plotHeatmap13.png')
		outfile14   = path.join(tmpdir, 'plotHeatmap14.png')
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
			plotHeatmap(m1, "%s", dendro = 'col', rows = c('ROW2', 'ROW3', 'ROW4', 'ROW5', 'ROW6'))
			plotHeatmap(m1, "%s", dendro = 'row', cols = c('COL2', 'COL3', 'COL4', 'COL5', 'COL6', 'COL7'))
			plotHeatmap(m1, "%s", dendro = FALSE, rows = c('ROW6', 'ROW7', 'ROW8', 'ROW9', 'ROW10'), cols = c('COL7', 'COL8', 'COL9', 'COL10'))
			""" % (outfile1, outfile2, outfile3, outfile4, outfile5, outfile6, outfile7, outfile8, outfile9, outfile10, outfile11, outfile12, outfile13, outfile14))
		rc = Popen(['Rscript', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		if rc != 0:
			self.fail('Script %s failed!' % scriptfile)

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
		rc = Popen(['Rscript', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
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
		rc = Popen(['Rscript', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
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
		rc = Popen(['Rscript', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)

	def testTxtSampleinfoPython(self):
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
			f.write("sample1\tp2\tgroup2\tbatch2\n")
			f.write("sample2\tp3\tgroup2\tbatch2\n")
			f.write("sample3\tp3\tgroup2\tbatch2\n")
		with open(scriptfile, 'w') as f:
			f.write(txt.sampleinfo.py)
			f.write("""
saminfo = txtSampleinfo("%s")
assert saminfo.rownames == ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6']
assert saminfo.colnames == ['Group']
assert saminfo.ncol == 1
assert saminfo.nrow == 6

saminfo = txtSampleinfo("%s")
assert saminfo.rownames == ['sample1', 'sample2', 'sample3', 'sample1', 'sample2', 'sample3']
assert saminfo.colnames == ['Patient', 'Group', 'Batch']
assert saminfo.ncol == 3
assert saminfo.nrow == 6

try:
	txtSampleinfo("%s")
	raise RuntimeError()
except ValueError:
	pass

""" % (groupfile2, groupfile3, groupfileE))

		rc = Popen(['python', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)

	def testTxtSampleinfoR(self):
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
				
				stopifnot( rownames(r2) != c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6') )
				stopifnot( r2$Group == c('group1', 'group1', 'group1', 'group2', 'group2', 'group2') )

				r3 = txtSampleinfo("%s")
				stopifnot( rownames(r3) != c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6') )
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

		rc = Popen(['Rscript', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)

	def testTxtFilter(self):
		infile     = path.join(tmpdir, 'txtfilter-in.txt')
		outfile1   = path.join(tmpdir, 'txtfilter-out1.txt')
		outfile2   = path.join(tmpdir, 'txtfilter-out2.txt')
		outfile3   = path.join(tmpdir, 'txtfilter-out3.txt')
		outfile4   = path.join(tmpdir, 'txtfilter-out4.txt')
		scriptfile = path.join(tmpdir, 'txtfilter.py')
		colnames   = ['a', 'b', 'c', 'd', 'e']
		rownames   = ['__r1', 'r2', 'r3']
		with open(infile, 'w') as f:
			f.write('\t'.join(colnames) + '\n')
			for rname in rownames:
				f.write(rname + '\t' + '\t'.join(random.sample(map(str, list(range(100))), 5)) + '\n')
		with open(scriptfile, 'w') as f:
			f.write(txt.filter.py)
			f.write('txtFilter("%s", "%s")\n' % (infile, outfile1))
			f.write('txtFilter("%s", "%s", [0,1,4], lambda x: True)\n' % (infile, outfile2))
			f.write('txtFilter("%s", "%s", [], lambda x: not x[0].startswith("__"))\n' % (infile, outfile3))
			f.write('txtFilter("%s", "%s", [], lambda x: True, False, 1)\n' % (infile, outfile4))

		rc = Popen(['python', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)

		def rcnames(outfile):
			with open(outfile) as f:
				rownames = []
				colnames = f.readline().strip().split('\t')
				for line in f:
					rownames.append(line.split('\t')[0])
			return colnames, rownames
		
		self.assertEqual(rcnames(outfile1), (colnames, rownames))
		self.assertEqual(rcnames(outfile2), (['a', 'd'], rownames))
		self.assertEqual(rcnames(outfile3), (colnames, ['r2', 'r3']))
		self.assertEqual(rcnames(outfile4)[1], ['r2', 'r3'])

	def testTxtTransformer(self):
		infile     = path.join(tmpdir, 'txttransformer-in.txt')
		outfile1   = path.join(tmpdir, 'txttransformer-out1.txt')
		scriptfile = path.join(tmpdir, 'txttransformer.py')
		colnames   = ['a', 'b', 'c', 'd', 'e']
		rownames   = ['__r1', 'r2', 'r3']
		with open(infile, 'w') as f:
			f.write('\t'.join(colnames) + '\n')
			for rname in rownames:
				f.write(rname + '\t' + '\t'.join(['1'] * 5) + '\n')
		with open(scriptfile, 'w') as f:
			f.write(txt.transform.py)
			f.write('txtTransform("%s", "%s", lambda parts: [p if i==0 else "a"+str(p) for i,p in enumerate(parts)], True)\n' % (infile, outfile1))

		rc = Popen(['python', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)
		with open(outfile1) as f:
			self.assertIn('a1', f.read())

	def testparams2CmdArgs(self):
		scriptfile = path.join(tmpdir, 'testparams2CmdArgs.py')
		with open(scriptfile, 'w') as f:
			f.write('from collections import OrderedDict\n')
			f.write(helpers.params2CmdArgs.py)
			f.write('params = OrderedDict([\n')
			f.write('	("a", True),\n')
			f.write('	("b", "c"),\n')
			f.write('	("-Xms128M -Xmx2G", True),\n')
			f.write('	("cde", 12)\n')
			f.write('])\n')
			f.write('print params2CmdArgs(params)\n')
			f.write('print params2CmdArgs(params, noq = ["cde"])\n')
			f.write('print params2CmdArgs(params, dash = "--", equal = "=", noq = ["cde"])\n')
			f.write('print params2CmdArgs(params, dash = "", equal = "=", noq = ["cde"])\n')
		stdout = check_output(['python', scriptfile], stderr=self.devnull)
		self.assertIn("-a -b c ---Xms128M -Xmx2G --cde=12", stdout.splitlines())
		self.assertIn("-a -b c ---Xms128M -Xmx2G --cde=12", stdout.splitlines())
		self.assertIn("a b=c -Xms128M -Xmx2G cde=12", stdout.splitlines())

	def testparams2CmdArgsR(self):
		scriptfile = path.join(tmpdir, 'testparams2CmdArgs.r')
		with open(scriptfile, 'w') as f:
			f.write(helpers.params2CmdArgs.r)
			f.write('params = list(\n')
			f.write('	a = TRUE,\n')
			f.write('	b = "c",\n')
			f.write('	cde = 12\n')
			f.write(')\n')
			f.write('params[["-Xms128M -Xmx2G"]] = T\n')
			f.write('write(params2CmdArgs(params), stdout())\n')
			f.write('write(params2CmdArgs(params, noq = c("cde")), stdout())\n')
			f.write('write(params2CmdArgs(params, dash = "--", equal = "=", noq = c("cde")), stdout())\n')
			f.write('write(params2CmdArgs(params, dash = "", equal = "=", noq = c("cde")), stdout())\n')
		stdout = check_output(['Rscript', scriptfile], stderr=self.devnull)
		self.assertIn("-a -b 'c' --cde='12' ---Xms128M -Xmx2G", stdout.splitlines())
		self.assertIn("-a -b 'c' --cde=12 ---Xms128M -Xmx2G", stdout.splitlines())
		self.assertIn("a b='c' cde=12 -Xms128M -Xmx2G", stdout.splitlines())

	def testMem2(self):
		scriptfile = path.join(tmpdir, 'testmem2.py')
		with open(scriptfile, 'w') as f:
			f.write(mem2.py)
			f.write('print mem2(38)\n')
			f.write('print mem2("20g", "M")\n')
			f.write('print mem2("2048000")\n')
			f.write('print mem2("2048000", "java")\n')
		stdout = check_output(['python', scriptfile], stderr=self.devnull)
		self.assertIn('38K', stdout.splitlines())
		self.assertIn('20480M', stdout.splitlines())
		self.assertIn('2000M', stdout.splitlines())
		self.assertIn('-Xms250M -Xmx2000M', stdout.splitlines())

	def testMem2R(self):
		scriptfile = path.join(tmpdir, 'testmem2.r')
		with open(scriptfile, 'w') as f:
			f.write(mem2.r)
			f.write('print(mem2(38))\n')
			f.write('print(mem2("20g", "M"))\n')
			f.write('print(mem2("2048000"))\n')
			f.write('print(mem2("2048000", "java"))\n')
		stdout = check_output(['Rscript', scriptfile], stderr=self.devnull)
		self.assertIn('[1] "38K"', stdout.splitlines())
		self.assertIn('[1] "20480M"', stdout.splitlines())
		self.assertIn('[1] "2000M"', stdout.splitlines())
		self.assertIn('[1] "-Xms250M -Xmx2000M"', stdout.splitlines())

	def testRuncmd(self):
		scriptfile = path.join(tmpdir, 'testRuncmd.py')
		with open(scriptfile, 'w') as f:
			f.write(runcmd.py)
			f.write('print runcmd("cmddosnotexists", quit=False)\n')
			f.write('print runcmd("python --version")\n')
		try:
			stdout = check_output(['python', scriptfile], stderr=self.devnull).splitlines()
		except Exception:
			stdout = ''
		self.assertEqual(['False', 'True'], stdout)
	
	def testRuncmdR(self):
		scriptfile = path.join(tmpdir, 'testRuncmd.r')
		with open(scriptfile, 'w') as f:
			f.write(runcmd.r)
			f.write('print(runcmd("python --version"))\n')
			f.write('print(runcmd("nosuchcmd", quit = F))\n')
		try:
			stdout = check_output(['Rscript', scriptfile], stderr=self.devnull).splitlines()
		except Exception:
			stdout = ''
		self.assertEqual(['[1] 0', '[1] 127'], stdout)

	def testPollingNon1st(self):
		scriptfile = path.join(tmpdir, 'testPollingNon1st.py')
		flagfile   = path.join(tmpdir, 'testPollingNon1st.py.flag')
		if path.exists(flagfile): remove(flagfile)
		if path.exists(flagfile+ '.error'): remove(flagfile + '.error')
		with open(scriptfile, 'w') as f:
			f.write(polling.non1st.py)
			f.write('from multiprocessing import Process, JoinableQueue as Queue\n')
			f.write('from time import time\n')
			f.write('def worker(q):\n')
			f.write('	while True:\n')
			f.write('		if q.empty(): break\n')
			f.write('		index, cmd = q.get()\n')
			f.write('		pollingNon1st(index, cmd, "%s", t = 1)\n' % flagfile)
			f.write('		print "Process:", index, "done!"\n')
			f.write('		q.task_done()\n')
			f.write('\n')
			f.write('sq = Queue()\n')
			f.write('for i in range(5):\n')
			f.write('	sq.put((i, "sleep 3"))\n')
			f.write('\n')
			f.write('t0 = time()\n')
			f.write('for i in range(5):\n')
			f.write('	t = Process(target = worker, args=(sq, ))\n')
			f.write('	t.daemon = True\n')
			f.write('	t.start ()\n')
			f.write('sq.join()\n')
			f.write('print int(time() - t0)\n')
		stdout = check_output(['python', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(stdout[-1], ['3', '4'])

	def testPollingNon1stR(self):
		scriptfile = path.join(tmpdir, 'testPollingNon1st.r')
		flagfile   = path.join(tmpdir, 'testPollingNon1st.r.flag')
		if path.exists(flagfile): remove(flagfile)
		if path.exists(flagfile+ '.error'): remove(flagfile + '.error')

		with open(scriptfile, 'w') as f:
			f.write(polling.non1st.r)
			f.write('library(doParallel)\n')
			f.write('t0 = proc.time()\n')
			f.write('cl <- makeCluster(5)\n')
			f.write('registerDoParallel(cl)\n')
			f.write('foreach(i=1:5) %dopar% {\n')
			f.write('	pollingNon1st(i-1, "sleep 4", "%s", t = 1)\n' % flagfile)
			f.write('}\n')
			f.write('stopCluster(cl) \n')
			f.write('cat(proc.time() - t0)\n')

		stdout = check_output(['Rscript', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(int(float(stdout[-1].split()[2])), [4, 5, 6])

	def testPollingAll(self):
		scriptfile = path.join(tmpdir, 'testPollingAll.py')
		for i in list(range(5)):
			odir = path.join(tmpdir, str(i), 'output')
			if not path.exists(odir): makedirs(odir)
			flag = path.join(odir, 'flag')
			eflg = path.join(odir, 'flag.error')
			if path.exists(flag): remove(flag)
			if path.exists(eflg): remove(eflg)

		with open(scriptfile, 'w') as f:
			f.write(polling.all.py)
			f.write('from multiprocessing import Process, JoinableQueue as Queue\n')
			f.write('from time import time\n')
			f.write('def worker(q):\n')
			f.write('	while True:\n')
			f.write('		if q.empty(): break\n')
			f.write('		index, cmd = q.get()\n')
			f.write('		pollingAll("%s", 5, index, cmd, "flag", t = 1)\n' % (tmpdir))
			f.write('		print "Process:", index, "done!"\n')
			f.write('		q.task_done()\n')
			f.write('\n')
			f.write('sq = Queue()\n')
			f.write('for i in range(5):\n')
			f.write('	sq.put((i, "sleep " + str(i)))\n')
			f.write('\n')
			f.write('t0 = time()\n')
			f.write('for i in range(5):\n')
			f.write('	t = Process(target = worker, args=(sq, ))\n')
			f.write('	t.daemon = True\n')
			f.write('	t.start ()\n')
			f.write('sq.join()\n')
			f.write('print int(time() - t0)\n')
		stdout = check_output(['python', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(stdout[-1], ['4', '5'])

	def testPollingAllR(self):
		scriptfile = path.join(tmpdir, 'testPollingAll.r')
		for i in list(range(5)):
			odir = path.join(tmpdir, str(i), 'output')
			if not path.exists(odir): makedirs(odir)
			flag = path.join(odir, 'flag')
			eflg = path.join(odir, 'flag.error')
			if path.exists(flag): remove(flag)
			if path.exists(eflg): remove(eflg)

		with open(scriptfile, 'w') as f:
			f.write(polling.all.r)
			f.write('library(doParallel)\n')
			f.write('t0 = proc.time()\n')
			f.write('cl <- makeCluster(5)\n')
			f.write('registerDoParallel(cl)\n')
			f.write('foreach(i=1:5) %dopar% {\n')
			f.write('	pollingAll("%s", 5, i-1, paste("sleep", i - 1), "flag", t = .1)\n' % tmpdir)
			f.write('}\n')
			f.write('stopCluster(cl) \n')
			f.write('cat(proc.time() - t0)\n')

		stdout = check_output(['Rscript', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(int(float(stdout[-1].split()[2])), [7, 5, 6])

	def testPollingFirst(self):
		scriptfile = path.join(tmpdir, 'testPollingFirst.py')
		for i in list(range(5)):
			odir = path.join(tmpdir, str(i), 'output')
			if not path.exists(odir): makedirs(odir)
			flag = path.join(odir, 'flag')
			eflg = path.join(odir, 'flag.error')
			if path.exists(flag): remove(flag)
			if path.exists(eflg): remove(eflg)

		with open(scriptfile, 'w') as f:
			f.write(polling.first.py)
			f.write('from multiprocessing import Process, JoinableQueue as Queue\n')
			f.write('from time import time\n')
			f.write('def worker(q):\n')
			f.write('	while True:\n')
			f.write('		if q.empty(): break\n')
			f.write('		index, cmd = q.get()\n')
			f.write('		pollingFirst("%s", 5, index, cmd, "flag", t = 1)\n' % (tmpdir))
			f.write('		print "Process:", index, "done!"\n')
			f.write('		q.task_done()\n')
			f.write('\n')
			f.write('sq = Queue()\n')
			f.write('for i in range(5):\n')
			f.write('	sq.put((i, "sleep " + str(i)))\n')
			f.write('\n')
			f.write('t0 = time()\n')
			f.write('for i in range(5):\n')
			f.write('	t = Process(target = worker, args=(sq, ))\n')
			f.write('	t.daemon = True\n')
			f.write('	t.start ()\n')
			f.write('sq.join()\n')
			f.write('print int(time() - t0)\n')
		stdout = check_output(['python', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(stdout[-1], ['4', '5'])

	def testPollingFirstR(self):
		scriptfile = path.join(tmpdir, 'testPollingFirst.r')
		for i in list(range(5)):
			odir = path.join(tmpdir, str(i), 'output')
			if not path.exists(odir): makedirs(odir)
			flag = path.join(odir, 'flag')
			eflg = path.join(odir, 'flag.error')
			if path.exists(flag): remove(flag)
			if path.exists(eflg): remove(eflg)

		with open(scriptfile, 'w') as f:
			f.write(polling.first.r)
			f.write('library(doParallel)\n')
			f.write('t0 = proc.time()\n')
			f.write('cl <- makeCluster(5)\n')
			f.write('registerDoParallel(cl)\n')
			f.write('foreach(i=1:5) %dopar% {\n')
			f.write('	pollingFirst("%s", 5, i-1, paste("sleep", i - 1), "flag", t = .1)\n' % tmpdir)
			f.write('}\n')
			f.write('stopCluster(cl) \n')
			f.write('cat(proc.time() - t0)\n')

		stdout = check_output(['Rscript', scriptfile], stderr=self.devnull).splitlines()
		self.assertIn(int(float(stdout[-1].split()[2])), [4, 5, 6])

	def testGeneNorm(self):
		scriptfile = path.join(tmpdir, 'genenorm.py')
		genefile   = path.join(tmpdir, 'gene.txt')
		retfile    = path.join(tmpdir, 'gene-ret.txt')
		with open(genefile, 'w') as f:
			f.write('CXCL1\n')
			f.write('EMMPRIN\n')
		with open(scriptfile, 'w') as f:
			f.write(genenorm.py)
			f.write("""
genenorm("%s", "%s")
			""" % (genefile, retfile))
		rc = Popen(['python', scriptfile], stdout=self.devnull, stderr=self.devnull).wait()
		self.assertEqual(rc, 0)


if __name__ == '__main__':
	unittest.main(failfast=True, verbosity = 2)