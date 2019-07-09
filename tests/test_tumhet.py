from pathlib import Path
import pytest
from remotedata import remotedata
from pyppl import PyPPL
from bioprocs.tumhet import pClonEvol, pPyClone
remotedata.manager.cachedir   = Path(__file__).parent / 'testdata'
remotedata.manager.conf.repos = 'pwwang/bioprocs-testdata'

def test_clonevol():
	pClonEvol1 = pClonEvol.copy()
	pClonEvol1.input = (
		remotedata.get('tumhet/aml1.txt'), remotedata.get('tumhet/aml1.sample.txt'))
	PyPPL().start(pClonEvol1).run()

def test_pyclone():
	vcfdir = remotedata.get('tumhet/SRR385940-individuals/')
	pPyClone1 = pPyClone.copy()
	pPyClone1.input = (
		list(vcfdir.glob('*.vcf.gz')),
		list(vcfdir.glob('*.vcf.gz')))
	PyPPL().start(pPyClone1).run()
