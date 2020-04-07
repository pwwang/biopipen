from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.stats import pSurvival, pCorr, pDiffCorr, pAdjust, pDeCov, pChow, pMediation
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
    assertInfile(pDiffCorr1.channel.outfile.get(), 'Compare2	V4	V8	0.941')
    assert Path(pDiffCorr1.channel.outdir.get() / 'Compare2-V4-V8.png').is_file()

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
    assertNotInfile(pDiffCorr2.channel.outfile.get(0), 'Compare1', 'Compare2')
    assertNotInfile(pDiffCorr2.channel.outfile.get(1), 'Compare1', 'Compare2')

def test_diffcorr_exhaust(rdata):
    pDiffCorr3 = pDiffCorr.copy()
    pDiffCorr3.input = (rdata.get('stats/corr.txt'),
                        rdata.get('stats/colfile.txt'))
    pDiffCorr3.args.nthread = 10
    PyPPL().start(pDiffCorr3).run()

def test_padj_header(rdata):
    pAdjust1 = pAdjust.copy()
    pAdjust1.input = [[rdata.get('stats/padj-header-in1.txt'),
                       rdata.get('stats/padj-header-in2.txt')]]
    pAdjust1.args.inopts.rnames = True
    pAdjust1.args.pcol = 1
    PyPPL().start(pAdjust1).run()
    assertInfile(pAdjust1.channel.outfile.get(), 'pval	Padj')
    assertInfile(pAdjust1.channel.outfile.get(), 'R1	1e-10	2e-10')
    assertInfile(pAdjust1.channel.outfile.get(), 'R4	1e-06	1.333')

def test_padj_noheader(rdata):
    pAdjust2 = pAdjust.copy()
    pAdjust2.input = [[rdata.get('stats/padj-noheader-in1.txt'),
                       rdata.get('stats/padj-noheader-in2.txt')]]
    pAdjust2.args.inopts.rnames = False
    pAdjust2.args.inopts.cnames = False
    pAdjust2.args.pcol = 2
    PyPPL().start(pAdjust2).run()
    assertNotInfile(pAdjust2.channel.outfile.get(), 'pval	Padj')
    assertInfile(pAdjust2.channel.outfile.get(), 'R1	1e-10	1e-10')
    assertInfile(pAdjust2.channel.outfile.get(), 'R2	1e-20	2e-20')

def test_decov(rdata):
    pDeCov1 = pDeCov.copy()
    pDeCov1.input = rdata.get('stats/corr.txt'), rdata.get('stats/cov.txt')
    PyPPL().start(pDeCov1).run()
    assertInfile(pDeCov1.channel.get(), 'V1	V2	V3')
    assertInfile(pDeCov1.channel.get(), 'S1	-14.37')
    assertInfile(pDeCov1.channel.get(), 'S2	12.48')
    assertInfile(pDeCov1.channel.get(), 'S10	6.48')

def test_chow(rdata):
    pChow1 = pChow.copy()
    pChow1.input = (rdata.get('stats/chow.txt'),
                    rdata.get('stats/chow-cols.txt'))
    pChow1.args.stacked = True
    pChow1.args.pval = 1.1
    PyPPL().start(pChow1).run()
    assertInfile(pChow1.channel.outfile.get(), 'P1:P2	Income	Savings	26')
    assertInfile(pChow1.channel.outfile.get(), '10.69') # Fstat

def test_mediation(rdata):
    pMediation1 = pMediation.copy()
    pMediation1.input = rdata.get('stats/corr.txt'), rdata.get('stats/medcases.txt')
    pMediation1.args.pval = 1
    PyPPL().start(pMediation1).run()
    assertInfile(pMediation1.channel.outfile.get(), 'Case1	S3	S2	S1	-0.077	-0.604	0.268')
    assertInfile(pMediation1.channel.outfile.get(), 'Case4	S8	S5	S9	-0.022	-0.293	0.856')
