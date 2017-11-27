import unittest
from os import path
from glob import glob
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config
from bioprocs.marray import pCeldir2Matrix, pBatchEffect, pMarrayDeg
from bioprocs.web import pDownloadGet
from bioprocs.matrix import pTxtFilter

class TestMarray (unittest.TestCase):

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
		pMakeAnno              = pTxtFilter.copy()
		pMakeAnno.depends      = pExtract
		pMakeAnno.input        = lambda ch: glob(path.join(ch.get(-1), '*.csv'))[0]
		pMakeAnno.args.skip    = 4
		pMakeAnno.args.header  = True
		pMakeAnno.args.cols    = [1,3]
		pMakeAnno.args.delimit = ','
		pMakeAnno.exdir        = getfile()
		PyPPL(config).start(pDownloadGet).run()
		
	def test1pExpdir2Matrix(self):
		pCeldir2Matrix.input         = [getfile()]
		pCeldir2Matrix.args.pattern  = "*.CEL"
		pCeldir2Matrix.args.cdffile  = getfile('miRNA-4_0-st-v1.cdf')
		pCeldir2Matrix.args.annofile = getfile('miRNA-4_0-st-v1.annotations.20160922.csv')
		pCeldir2Matrix.args.boxplot  = True
		pCeldir2Matrix.args.heatmap  = True
		pCeldir2Matrix.args.histplot = True
		pCeldir2Matrix.args.norm     = 'rma'

		PyPPL(config).start(pCeldir2Matrix).run()
		procOK(pCeldir2Matrix, 'celdir2mat.txt', self)
	
	def test2pBatchEffect(self):
		pBatchEffect.input = (getfile('celdir2mat.txt', input = False), getfile('saminfo.txt'))
		pBatchEffect.args.boxplot  = True
		pBatchEffect.args.heatmap  = True
		pBatchEffect.args.histplot = True
		PyPPL(config).start(pBatchEffect).run()
		procOK(pBatchEffect, 'batcheffect.txt', self)

	def test3pDeg(self):

		pDegPaired        = pMarrayDeg.copy()
		pMarrayDeg.input        = (getfile('celdir2mat.txt', input = False), getfile('saminfo.txt'))
		pMarrayDeg.args.heatmap = True
		pMarrayDeg.args.maplot  = True

		pDegPaired.args.heatmap = True
		pDegPaired.args.paired  = True
		pDegPaired.input        = (getfile('celdir2mat.txt', input = False), getfile('saminfo-paired.txt'))
		PyPPL(config).start(pMarrayDeg, pDegPaired).run()
		procOK(pMarrayDeg, 'deg.txt', self)
		procOK(pDegPaired, 'deg-paired.txt', self)

if __name__ == '__main__':
	unittest.main(failfast=True)