"""RNA-seq data analysis"""

from ..core.proc import Proc
from ..core.config import config

class UnitConversion(Proc):
    """Convert expression value units back and forth"""
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.rscript
    envs = {
        "infmt": "matrix", # or rds
        "inunit": None,
        "outunit": None,
        "refexon": config.ref.refexon,
        "meanfl": None,
        "inlog2p": False,
        "outlog2p": False,
    }
    script = "file://../scripts/rnaseq/UnitConversion.R"
