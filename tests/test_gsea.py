from pathlib import Path
import pytest
from pyppl import PyPPL
from biopipen.gsea import pEnrichr

def test_enrichr():
	pEnrichr1 = pEnrichr.copy()
