from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.stats import pSurvival
from . import remotedata

def test_survival():
	pSurvival1 = pSurvival.copy()
	pSurvival1.input = [remotedata.get('stats/survival_cat.txt')]
	PyPPL().start(pSurvival1).run()
