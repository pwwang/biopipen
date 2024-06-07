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
        tabix: Path to tabix, used to index infile/outfile
        ncores (type=int): Number of cores (`--threads`) to use
        keep: Whether we should keep the filtered variants or not.
            If True, the filtered variants will be kept in the output file, but
            with a new FILTER.
        includes: and
        excludes: include/exclude only sites for which EXPRESSION is true.
            See: <https://samtools.github.io/bcftools/bcftools.html#expressions>
            If provided, `envs.include/exclude` will be ignored.
            If `str`/`list` used, The filter names will be `Filter_<type>_<index>`.
            A dict is used where keys are filter names and values are expressions
        gz (flag): Whether to gzip the output file
        index (flag): Whether to index the output file (tbi) (`envs.gz` forced to True)
        <more>: Other arguments for `bcftools filter`
            See also <https://samtools.github.io/bcftools/bcftools.html#filter>
    """
    input = "infile:file"
    output = (
        "outfile:file:{{in.infile | stem: 'gz'}}.vcf"
        "{{'.gz' if envs.index or envs.gz else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
        "keep": True,
        "includes": None,
        "excludes": None,
        "gz": True,
        "index": True,
    }
    script = "file://../scripts/bcftools/BcftoolsFilter.py"


class BcftoolsSort(Proc):
    """Sort VCF files using `bcftools sort`.

    `bcftools sort` is used to sort VCF files by chromosome and position based on the
    order of contigs in the header.

    Here we provide a chrsize file to first sort the contigs in the header and then
    sort the VCF file using `bcftools sort`.

    Input:
        infile: The input VCF file

    Output:
        outfile: The sorted VCF file.

    Envs:
        bcftools: Path to bcftools
        tabix: Path to tabix, used to index infile/outfile
        ncores (type=int): Number of cores (`--threads`) to use
        gz (flag): Whether to gzip the output file
        index (flag): Whether to index the output file (tbi) (`envs.gz` forced to True)
        chrsize: The chromosome size file, from which the chromosome order is used
            to sort the contig in the header first.
            If not provided, `bcftools sort` will be used directly.
        notfound (choice): What if the contig in the VCF file is not found in the
            `chrsize` file.
            - error: Report error
            - remove: Remove the contig from the header.
                Note that if there are records with the removed contig, an error will
                be raised by `bcftools sort`
            - start: Move the contig to the start of the contigs from `chrsize`
            - end: Move the contig to the end of the contigs from `chrsize`
        <more>: Other arguments for `bcftools sort`. For example `max_mem`.
            See also <https://samtools.github.io/bcftools/bcftools.html#sort>
    """
    input = "infile:file"
    output = (
        "outfile:file:{{in.infile | stem: 'gz'}}.vcf"
        "{{'.gz' if envs.index or envs.gz else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
        "chrsize": config.ref.chrsize,
        "notfound": "remove",
        "gz": True,
        "index": True,
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
