from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.gsea import pEnrichr
from . import remotedata

def test_enrichr():
	pEnrichr1 = pEnrichr.copy()
