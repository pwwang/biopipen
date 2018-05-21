import unittest, helpers, testly
from os import path
from glob import glob
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config
from bioprocs.marray import pCELdir2Matrix, pBatchEffect, pMArrayDEG
from bioprocs.web import pDownloadGet
from bioprocs.tsv import pTsv

class TestMarray (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestMarray')

	@classmethod
	def setUpClass(self):
		pDownloadGet.input = [
			'https://www.dropbox.com/s/gx6oi0zx72dolgu/genechip_4-0_array_sample_data.zip?dl=1', 'https://www.dropbox.com/s/vlab9fzfskg8pve/miRNA-4_0-st-v1_CDF.zip?dl=1', 'https://www.dropbox.com/s/4wtjyb38tptpxuu/miRNA-4_0-st-v1-annotations-20160922-csv.zip?dl=1'
		]
		pDownloadGet.forks = 3

		pExtract         = Proc(desc = 'Extract files')
		pExtract.depends = pDownloadGet
		pExtract.forks   = 3
		pExtract.input   = "infile:file"
		pExtract.output  = "outdir:dir:celdata"
		pExtract.script  = """
		unzip {{in.infile}} -d {{out.outdir}}
		if ls {{out.outdir}}/*/* 2>/dev/null; then
			mv {{out.outdir}}/*/* {{out.outdir}}
		fi
		"""
		pExtract.exdir         = getfile()
		pExtract.expart        = ["{{out.outdir}}/*.CEL", "{{out.outdir}}/*.cdf"]
		pMakeAnno              = pTsv.copy()
		pMakeAnno.depends      = pExtract
		pMakeAnno.input        = lambda ch: glob(path.join(ch.get(-1), '*.csv'))[0]
		pMakeAnno.args.inopts.skip = 4
		pMakeAnno.args.inopts.ftype = 'head'
		pMakeAnno.args.inopts.delimit = '","'
		pMakeAnno.args.outopts.head  = True
		pMakeAnno.args.outopts.cnames = ['Probe Set Name', 'Transcript ID(Array Design)']
		pMakeAnno.exdir        = self.indir
		PyPPL(config).start(pDownloadGet).run()

	def test1pExpdir2Matrix(self):
		pCELdir2Matrix.input               = [self.indir]
		pCELdir2Matrix.args.pattern        = "*.CEL"
		pCELdir2Matrix.args.cdffile        = path.join(self.indir, 'miRNA-4_0-st-v1.cdf')
		pCELdir2Matrix.args.annofile       = path.join(self.indir, 'miRNA-4_0-st-v1.annotations.20160922.tsv')
		pCELdir2Matrix.args.plot.boxplot   = True
		pCELdir2Matrix.args.plot.heatmap   = True
		pCELdir2Matrix.args.plot.histogram = True
		pCELdir2Matrix.args.fn2sample      = 'function(fn) unlist(strsplit(fn, "_(", fixed=T))[1]'
		pCELdir2Matrix.args.norm           = 'rma'

		PyPPL(config).start(pCELdir2Matrix).run()
		self.assertFileEqual(pCELdir2Matrix.channel.get(0), path.join(self.outdir, 'celdir2mat.txt'))

	def test2pBatchEffect(self):
		pBatchEffect.input = (getfile('celdir2mat.txt', input = False), getfile('saminfo.txt'))
		pBatchEffect.args.plot.boxplot  = True
		pBatchEffect.args.plot.heatmap  = True
		pBatchEffect.args.plot.histogram = True
		PyPPL(config).start(pBatchEffect).run()
		self.assertFileEqual(pBatchEffect.channel.get(0), path.join(self.outdir, 'batcheffect.txt'))

	def dataProvider_test3pDeg(self):
		efile1   = path.join(self.outdir, 'celdir2mat.txt')
		gfile1   = path.join(self.indir, 'saminfo.txt')
		outfile1 = path.join(self.outdir, 'deg.txt')
		efile2   = path.join(self.outdir, 'celdir2mat.txt')
		gfile2   = path.join(self.indir, 'saminfo-paired.txt')
		outfile2 = path.join(self.outdir, 'deg-paired.txt')
		yield 't1', efile1, gfile1, outfile1
		yield 't2', efile1, gfile2, outfile2

	def test3pDeg(self, tag, efile, gfile, outfile):
		pMArrayDEGTest = pMArrayDEG.copy(tag = tag)
		pMArrayDEGTest.input = (efile, gfile)
		pMArrayDEGTest.args.plot.heatmap = True
		pMArrayDEGTest.args.plot.maplot  = True
		pMArrayDEGTest.args.plot.mdsplot = True
		pMArrayDEGTest.args.plot.volplot = True
		PyPPL(config).start(pMArrayDEGTest).run()
		self.assertFileEqual(pMArrayDEGTest.channel.get(), outfile)


if __name__ == '__main__':
	testly.main(failfast=True)
