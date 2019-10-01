from pathlib import Path
import pytest
from pyppl import PyPPL, Box
from bioprocs.tumhet import pClonEvol, pPyClone, pPyClone2ClonEvol, pAllFIT
from . import remotedata

def test_clonevol():
	pClonEvol1 = pClonEvol.copy()
	pClonEvol1.input = (
		remotedata.get('tumhet/aml1.txt'), remotedata.get('tumhet/aml1.sample.txt'))
	PyPPL().start(pClonEvol1).run()

def test_pyclone():
	vcfdir = remotedata.get('tumhet/SRR385940-individuals/')
	pPyClone1 = pPyClone.copy()
	pPyClone2ClonEvol1 = pPyClone2ClonEvol.copy()
	pClonEvol2 = pClonEvol.copy()

	pPyClone1.input = (
		list(vcfdir.glob('*.vcf.gz')),
		list(vcfdir.glob('*.vcf.gz')))
	pPyClone2ClonEvol1.depends = pPyClone1
	pClonEvol2.depends = pPyClone2ClonEvol1
	pClonEvol2.input = lambda ch: ch.cbind(remotedata.get('tumhet/SRR385940.sample.txt'))
	PyPPL().start(pPyClone1).run().report()

def test_allfit():
	from remotedata import GithubRemoteData
	remotedata2 = GithubRemoteData(Box(
		source = 'github',
		cachedir = Path(__file__).parent / 'testdata',
		repos = 'KhiabanianLab/All-FIT'
	))
	infile = remotedata2.get('test/input/sampleFile1.xls')
	pAllFIT1 = pAllFIT.copy()
	pAllFIT1.input = [infile]
	PyPPL().start(pAllFIT1).run()

