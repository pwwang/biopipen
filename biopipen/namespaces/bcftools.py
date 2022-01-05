"""handling VCF files using bcftools"""
from ..core.proc import Proc
from ..core.config import config

class BcftoolsAnnotate(Proc):
    """Add or remove annotations from VCF files

    Input:
        infile: The input VCF file
        annfile: The annotation file

    Output:
        outfile: The annotated VCF file

    Envs:
        bcftools: Path to bcftools
        tabix: Path to tabix, used to index infile and annfile
        annfile: The annotation file. If `in.annfile` is provided,
            this is ignored
        ncores: Number of cores (`--nthread`) to use
        cols: Overwrite `-c/--columns`
        header: Headers to be added
        args: Other arguments for `bcftools annotate`
    """
    input = "infile:file, annfile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "annfile": "",
        "cols": [],
        "header": [],
        "args": {},
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/bcftools/BcftoolsAnnotate.py"
