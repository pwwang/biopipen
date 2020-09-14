from pathlib import Path
import pytest
from diot import Diot
from pyppl import PyPPL
from biopipen.tumhet import pClonEvol, pPyClone, pPyClone2ClonEvol, pAllFIT, pTheta

def test_clonevol(rdata):
    pClonEvol1 = pClonEvol.copy()
    pClonEvol1.input = (
        rdata.get('tumhet/aml1.txt'), rdata.get('tumhet/aml1.sample.txt'))
    PyPPL().start(pClonEvol1).run()

def test_pyclone(rdata):
    vcf1 = str(rdata.get('tumhet/SRR385940-individuals/SRR385940-NA12156.vcf.gz'))
    vcf2 = str(rdata.get('tumhet/SRR385940-individuals/SRR385940-NA12878.vcf.gz'))
    vcf3 = str(rdata.get('tumhet/SRR385940-individuals/SRR385940-NA18507.vcf.gz'))
    vcf4 = str(rdata.get('tumhet/SRR385940-individuals/SRR385940-NA19240.vcf.gz'))
    pPyClone1 = pPyClone.copy()
    pPyClone2ClonEvol1 = pPyClone2ClonEvol.copy()
    pClonEvol2 = pClonEvol.copy()

    pPyClone1.input = (','.join([vcf1, vcf2, vcf3, vcf4]), ) * 2
    pPyClone2ClonEvol1.depends = pPyClone1
    pClonEvol2.depends = pPyClone2ClonEvol1
    pClonEvol2.input = lambda ch: ch.cbind(rdata.get('tumhet/SRR385940.sample.txt'))
    PyPPL().start(pPyClone1).run().report()

def test_allfit(rdata):
    from remotedata import GithubRemoteData
    rdata2 = GithubRemoteData(Diot(
        source = 'github',
        cachedir = Path(__file__).parent / 'testdata',
        repos = 'KhiabanianLab/All-FIT'
    ))
    infile = rdata2.get('test/input/sampleFile1.xls')
    pAllFIT1 = pAllFIT.copy()
    pAllFIT1.input = [infile]
    PyPPL().start(pAllFIT1).run()


def test_theta(config):
    from remotedata import remotedata
    rd3 = remotedata(Diot(
        source   = config.source,
        cachedir = config.cachedir,
        repos    = 'raphael-group/THetA'
    ))
    pTheta1 = pTheta.copy()
    pTheta1.input = (
        rd3.get('example/Example.intervals'),
        rd3.get('example/TUMOR_SNP.formatted.txt'),
        rd3.get('example/NORMAL_SNP.formatted.txt'),
    )
    PyPPL().start(pTheta1).run()
