from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.stats import pSurvival, pCorr
from . import assertInfile

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
