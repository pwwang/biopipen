from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.stats import pSurvival, pCorr, pDiffCorr
from . import assertInfile, assertNotInfile

def test_survival(rdata):
    pSurvival1 = pSurvival.copy()
    pSurvival1.input = [rdata.get('stats/survival_cat.txt')]
    PyPPL().start(pSurvival1).run()

def test_corr(rdata):
    pCorr1           = pCorr.copy()
    pCorr1.input     = [rdata.get('stats/corr.txt')]
    pCorr1.args.pval = True
    pCorr1.args.plot = True
    PyPPL().start(pCorr1).run()
    assertInfile(str(pCorr1.channel.outfile.get())[:-3] + 'pairs.txt', 'V1	V1	1')

def test_diffcorr(rdata):
    pDiffCorr1 = pDiffCorr.copy()
    pDiffCorr1.input = (rdata.get('stats/corr.txt'),
                        rdata.get('stats/colfile.txt'),
                        rdata.get('stats/case.txt'),
                        rdata.get('stats/rowfile.txt'))
    pDiffCorr1.args.stacked = False
    pDiffCorr1.args.plot = True
    pDiffCorr1.args.nthread = 10
    PyPPL().start(pDiffCorr1).run()
    assertInfile(pDiffCorr1.channel.outfile.get(), 'Case1	V8	V4	0.609')
    assert Path(pDiffCorr1.channel.outdir.get() / 'Case1-V8-V4.png').is_file()

def test_diffcorr_stacked(rdata):
    pDiffCorr2 = pDiffCorr.copy()
    pDiffCorr2.input = [(rdata.get('stats/corr.txt'),
                         rdata.get('stats/colfile_stacked.txt'),
                         rdata.get('stats/case_stacked.txt'),
                         rdata.get('stats/rowfile_stacked.txt')),
                        (rdata.get('stats/corr.txt'),
                         rdata.get('stats/colfile_stacked.txt'),
                         rdata.get('stats/case_stacked_with_row_groups.txt'),
                         rdata.get('stats/rowfile_stacked.txt'))]
    pDiffCorr2.args.stacked = True
    pDiffCorr2.args.nthread = 10
    PyPPL().start(pDiffCorr2).run()
    assertNotInfile(pDiffCorr2.channel.outfile.get(0), 'Case1', 'Case2')
    assertNotInfile(pDiffCorr2.channel.outfile.get(1), 'Case1', 'Case2')

def test_diffcorr_exhaust(rdata):
    pDiffCorr3 = pDiffCorr.copy()
    pDiffCorr3.input = (rdata.get('stats/corr.txt'),
                        rdata.get('stats/colfile.txt'))
    pDiffCorr3.args.nthread = 10
    PyPPL().start(pDiffCorr3).run()