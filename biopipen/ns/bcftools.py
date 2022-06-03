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


class BcftoolsFilter(Proc):
    """Apply fixed threshold filters to VCF files

    Input:
        infile: The input VCF file

    Output:
        outfile: The filtered VCF file. If the `in.infile` is gzipped, this is
            gzipped as well.

    Envs:
        bcftools: Path to bcftools
        ncores: Number of cores (`--nthread`) to use
        keep: Whether we should keep the filtered variants or not.
        args: Other arguments for `bcftools annotate`
        ncores: `nthread`
        tmpdir: Path to save the intermediate files
            Since the filters need to be applied one by one by bcftools
        includes: and
        excludes: include/exclude only sites for which EXPRESSION is true.
            - See: https://samtools.github.io/bcftools/bcftools.html#expressions
            - If provided, `envs.args.include/exclude` will be ignored.
            - If `str`/`list` used, The filter names will be `Filter%d`
            - A dict is used when keys are filter names and values are
              expressions
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "keep": True,
        "ncores": config.misc.ncores,
        "includes": None,
        "excludes": None,
        "tmpdir": config.path.tmpdir,
        "args": {},
    }
    script = "file://../scripts/bcftools/BcftoolsFilter.py"
