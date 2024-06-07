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
            See: https://samtools.github.io/bcftools/bcftools.html#expressions
            If provided, `envs.args.include/exclude` will be ignored.
            If `str`/`list` used, The filter names will be `Filter%d`
            A dict is used when keys are filter names and values are expressions
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


class BcftoolsSort(Proc):
    """Sort VCF files

    Input:
        infile: The input VCF file

    Output:
        outfile: The sorted VCF file.

    Envs:
        bcftools: Path to bcftools
        gz: Whether to gzip the output file
        index: Whether to index the output file (tbi) (`envs.gz` forced to True)
        tmpdir: Path to save the intermediate files
        args: Other arguments for `bcftools sort`. For example `max-mem`.
            See also https://samtools.github.io/bcftools/bcftools.html#sort
    """
    input = "infile:file"
    output = (
        "outfile:file:{{in.infile | stem0}}.vcf"
        "{% if envs.gz or envs.index %}.gz{% endif %}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "gz": True,
        "index": True,
        "tmpdir": config.path.tmpdir,
        "args": {},
    }
    script = "file://../scripts/bcftools/BcftoolsSort.py"


class BcftoolsView(Proc):
    """View, subset and filter VCF files by position and filtering expression.

    Also convert between VCF and BCF.

    Input:
        infile: The input VCF file
        regions_file: The region file used to subset the input VCF file.
        samples_file: The samples file used to subset the input VCF file.

    Output:
        outfile: The output VCF file.

    Envs:
        bcftools: Path to bcftools
        tabix: Path to tabix, used to index infile/outfile
        ncores (type=int): Number of cores (`--threads`) to use
        regions_file: The region file used to subset the input VCF file.
            If `in.regions_file` is provided, this is ignored.
        samples_file: The samples file used to subset the input VCF file.
            If `in.samples_file` is provided, this is ignored.
        gz (flag): Whether to gzip the output file
        index (flag): Whether to index the output file (tbi) (`envs.gz` forced to True)
        <more>: Other arguments for `bcftools view`.
            See also https://samtools.github.io/bcftools/bcftools.html#view
            Note that the underscore `_` will be replaced with dash `-` in the
            argument name.
    """
    input = "infile:file, regions_file:file, samples_file:file"
    output = (
        "outfile:file:{{in.infile | stem: 'gz'}}.vcf"
        "{{'.gz' if envs.index or envs.gz else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
        "regions_file": None,
        "samples_file": None,
        "gz": True,
        "index": True,
    }
    script = "file://../scripts/bcftools/BcftoolsView.py"
