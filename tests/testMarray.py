import unittest, addPath

import random
from os import makedirs, path
from tempfile import gettempdir
from pyppl import PyPPL, Channel, Box, Proc
from bioprocs.marray import pCeldir2Matrix, pBatchEffect, pMarrayDeg
from bioprocs.web import pDownloadGet


unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestRnaseq (unittest.TestCase):

	wdir = path.realpath(path.join(path.dirname(path.dirname(__file__)), 'workdir'))
	data = {}

	@unittest.skipIf(path.exists(path.join(path.dirname(path.dirname(__file__)), 'workdir', 'celdata')), 'Data already generated.')
	def test0genData(self):
		pDownloadGet.input = ['https://www.dropbox.com/s/gx6oi0zx72dolgu/genechip_4-0_array_sample_data.zip?dl=1', 'https://www.dropbox.com/s/vlab9fzfskg8pve/miRNA-4_0-st-v1_CDF.zip?dl=1', 'https://www.dropbox.com/s/4wtjyb38tptpxuu/miRNA-4_0-st-v1-annotations-20160922-csv.zip?dl=1']
		pDownloadGet.forks   = 3

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
		pExtract.exdir   = path.join(self.wdir, 'celdata')
		pExtract.expart  = ["{{out.outdir}}/*.CEL", "{{out.outdir}}/*.cdf", "{{out.outdir}}/*.csv"]

		PyPPL().start(pDownloadGet).run()


	def test1pExpdir2Matrix(self):
		pCeldir2Matrix.input         = [path.join(self.wdir, 'celdata')]
		pCeldir2Matrix.args.pattern  = "*.CEL"
		pCeldir2Matrix.args.cdffile  = path.join(self.wdir, 'celdata', 'miRNA-4_0-st-v1.cdf')
		pCeldir2Matrix.args.annofile = path.join(self.wdir, 'celdata', 'miRNA-4_0-st-v1.annotations.20160922.csv')
		pCeldir2Matrix.args.boxplot  = True
		pCeldir2Matrix.args.heatmap  = True
		pCeldir2Matrix.args.histplot = True
		pCeldir2Matrix.args.norm     = 'rma'

		PyPPL().start(pCeldir2Matrix).run()
		self.data['expfile'] = pCeldir2Matrix.channel.outfile.flatten()[0]

	def test2pBatchEffect(self):
		pBatch = Proc(desc = 'Get batch')
		pBatch.input = {"infile:file": self.data['expfile']}
		pBatch.output = "outfile:file:group.txt"
		pBatch.script = """
		echo -e "Group\\tBatch" > {{out.outfile}}
		s=$(head -1 {{in.infile}} | tr "\\t" "\\n")
		echo "$s" | grep Brain | while read line; do
			echo -e "$line\\tBrain\\tBatch1" >> {{out.outfile}}
		done
		echo "$s" | grep Lung | while read line; do
			echo -e "$line\\tLung\\tBatch2" >> {{out.outfile}}
		done
		"""

		pBatchEffect.depends = pBatch
		pBatchEffect.input = lambda ch: ch.insert(0, self.data['expfile'])
		pBatchEffect.args.boxplot  = True
		pBatchEffect.args.heatmap  = True
		pBatchEffect.args.histplot = True
		PyPPL().start(pBatch).run()
		self.data['gfile'] = pBatch.channel.outfile.flatten()[0]
		self.data['expfile'] = pBatchEffect.channel.outfile.flatten()[0]

	def test3pDeg(self):

		pDegPaired        = pMarrayDeg.copy()
		pMarrayDeg.input        = (self.data['expfile'], self.data['gfile'])
		pMarrayDeg.args.heatmap = True
		pMarrayDeg.args.maplot  = True

		pDegPaired.args.heatmap = True
		pDegPaired.args.paired  = True
		pDegPaired.input        = (self.data['expfile'], self.data['gfile'])
		PyPPL().start(pMarrayDeg, pDegPaired).run()
		#PyPPL().start(pDeg).run()


if __name__ == '__main__':
	unittest.main(failfast=True)