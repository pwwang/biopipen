"""Tools to handle VCF files"""
from ..core.proc import Proc
from ..core.config import config


class VcfLiftOver(Proc):
    """Liftover a VCF file using GATK

    Input:
        invcf: The input VCF file

    Output:
        outvcf: The output VCF file

    Envs:
        gatk: The path to gatk4, which should be installed via conda
        chain: The map chain file for liftover
        tmpdir: Directory for temporary storage of working files
        args: Other CLI arguments for `gatk LiftoverVcf`
    """
    input = "invcf:file"
    output = "outvcf:file:{{in.invcf | basename}}"
    envs = {
        "gatk": config.exe.gatk4,
        "chain": config.path.liftover_chain,
        "tmpdir": config.path.tmpdir,
        "reffa": config.ref.reffa,
        "args": {},
    }
    lang = config.lang.bash
    script = "file://../scripts/vcf/VcfLiftOver.sh"


class VcfFilter(Proc):
    """Filter records in vcf file

    Input:
        invcf: The input vcf file, could be bgzipped.

    Output:
        outfile: The filtered vcf file. If `in.invcf` is bgzipped, then
            this will be bgzipped.

    Envs:
        filters: A dict of filters with keys the filter names.
            >>> # Typically
            >>> lambda variant: <expression>
            Things to notice
            1. Filters should return `False` to get variant filtered out
            2. See https://brentp.github.io/cyvcf2/docstrings.html#cyvcf2.cyvcf2.Variant
            For what you can do with the variant
            3. The filter python functions should be in string representation
            4. Builtin filters can have parameters `{"QUAL": 30}`
            5. List of builtin filters. Specify them like: `{"FILTER": params}`
            `SNPONLY`: keeps only SNPs (`{"SNPONLY": False}` to filter SNPs out)
            `QUAL`: keeps variants with QUAL>=param (`{"QUAL": (30, False)}`)
            to keep only variants with QUAL<30
        filter_descs: Descriptions for the filters. Will be saved to the header
            of the output vcf file
        helper: Some helper code for the filters
        keep: Keep the variants not passing the filters?
    """  # noqa: E501
    input = "invcf:file"
    output = "outfile:file:{{in.invcf | basename}}"
    lang = config.lang.python
    envs = {
        "filters": {},
        "keep": True,
        "helper": "",
        "filter_descs": {},
    }
    script = "file://../scripts/vcf/VcfFilter.py"


class VcfIndex(Proc):
    """Index VCF files. If they are already index, use the index files

    Input:
        infile: The input VCF file

    Output:
        outfile: The output VCF file (bgzipped)
        outidx: The index file of the output VCF file

    Envs:
        tabix: Path to tabix
    """
    input = "infile:file"
    output = """
        {%- if in.infile.endswith(".gz") %}
            outfile:file:{{in.infile | basename}},
            outidx:file:{{in.infile | basename | append: ".tbi"}}
        {%- else -%}
            outfile:file:{{in.infile | basename | append: ".gz"}},
            outidx:file:{{in.infile | basename | append: ".gz.tbi"}}
        {% endif -%}
    """
    lang = config.lang.python
    envs = {
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/vcf/VcfIndex.py"


class Vcf2Bed(Proc):
    """Convert Vcf file to Bed file

    Input:
        infile: The vcf file

    Output:
        outfile: The converted bed file

    Envs:
        inbase: The coordinate base of the vcf file
        outbase: The coordinate base of the base file

    Requires:
        cyvcf2:
            - check: {{proc.lang}} -c "import cyvcf2"
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem0}}.bed"
    lang = config.lang.python
    envs = {"inbase": 1, "outbase": 0}
    script = "file://../scripts/vcf/Vcf2Bed.py"


class VcfDownSample(Proc):
    """Down-sample VCF files to keep only a subset of variants in there

    Input:
        infile: The input VCF file

    Output:
        outfile: The output VCF file with subet variants
            Gzipped if `in.infile` is gzipped

    Envs:
        n: Fraction/Number of variants to keep
            If `n > 1`, it is the number.
            If `n <= 1`, it is the fraction.
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    envs = {"n": 0}
    lang = config.lang.bash
    script = "file://../scripts/vcf/VcfDownSample.sh"


class VcfSplitSamples(Proc):
    """Split a VCF file into multiple VCF files, one for each sample

    Input:
        infile: The input VCF file

    Output:
        outdir: The output directory containing the split VCF files

    Envs:
        bcftools: Path to bcftools
        gz: Gzip the output VCF files? Has to be True if `envs.index` is True
        index: Index the output VCF files?
        ncores: Number of cores, used to extract samples, but not to index
        private: Keep sites where only the sample carries an non-ref allele.
            That means, sites with genotypes like `0/0` will be removed.
    """
    input = "infile:file"
    output = "outdir:dir:{{in.infile | stem}}.splitsamples"
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "gz": True,
        "index": True,
        "ncores": config.misc.ncores,
        "private": True,
    }
    script = "file://../scripts/vcf/VcfSplitSamples.py"


class VcfIntersect(Proc):
    """Find variants in both VCF files

    Input:
        infile1: The first VCF file
        infile2: The second VCF file

    Output:
        outfile: The output VCF file with subet variants in both files

    Envs:
        bcftools: Path to bcftools
        gz: Gzip the output VCF files? Has to be True if `envs.index` is True
        index: Index the output VCF files?
        keep_as: Keep the variants as presented in the first (0) or
            the second (1) file?
        collapse: How to match the variants in the two files? Will be passed to
            `bcftools isec -c` option. See also
            https://samtools.github.io/bcftools/bcftools.html#common_options
            - none: only records with identical REF and ALT alleles are
                compatible
            - some: only records where some subset of ALT alleles match are
                compatible
            - all: all records are compatible, regardless of whether the ALT
                alleles match or not.
            - snps: any SNP records are compatible, regardless of whether the
                ALT alleles match or not.
            - indels: any indel records are compatible, regardless of whether
                the ALT alleles match or not.
            - both: abbreviates `snps` and `indels`
            - id: only records with identical ID are compatible
    """
    input = "infile1:file, infile2:file"
    output = """
        outfile:file:{{in.infile1 | stem0}}.intersect.{{in.infile2 | stem0}}.vcf
        {%- if envs.gz -%}.gz{%- endif -%}
    """
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "gz": True,
        "index": True,
        "collapse": "all",
    }
    script = "file://../scripts/vcf/VcfIntersect.py"


class VcfFix(Proc):
    """Fix some issues with VCF files

    Input:
        infile: The input VCF file

    Output:
        outfile: The output VCF file

    Envs:
        fixes: A list of fixes to apply.
            Each one is a dict with keys `kind`, `id`, `regex` and `fix`
            `kind`: The kind of fix. Including
                `filter` the FILTERs in the header,
                `info` the info INFOs in the header,
                `contig` the contig lines in the header
                `format` the FORMATs in the header,
                `colnames` the column names in the header
                `header` general header item
                `variant` the variants
                `None` matches everything
            `id`: The ID the match. If `kind` is `filter`, `info`, `contig`
                or `format`, then it matches the `ID` of the item. If `kind`
                is `variant`, then it matches the `ID` of the variant.
                If a list is given, then it matches any of the IDs in the list.
            `regex`: The regular expression to match. When `id` is given,
                this is ignored.
            `append`: Whether to append a record instead of to replace an
                existing one. When it is True, `kind` has to not be `None`
            `fix`: The fix to apply in the format of a lambda function
                (in string), with a single argument.
                The function should either return a string (raw representation)
                for the record, the record itself, `None`, `False`.
                If `None` is returned, the original record is used. if `False`,
                the record is removed.
                If `append` is `True`, then the function should either return
                a string or an object. And the argument is `None`
                The argument is a different object based on different `kind`s.
                When `kind` is `None`, the argument is the plain line of the
                record with line ending.
                When `kind` is `info` or `format`, the record is a dict with
                keys `ID`, `Description`, `Type` and `Number`.
                When `kind` is `filter`, the record is a dict with keys
                `ID` and `Description`.
                When `kind` is `contig`, the record is a dict with keys `ID`
                and `length`.
                When `kind` is `header`, the record is a dict
                with `key` the name of the header and `value` the value of the
                header.
                When `kind` is `colnames`, the record is a list of column names.
                When `kind` is `variant`, the record is a dict with
                keys `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`,
                `FORMAT` and `SAMPLES`. `INFO` is a dict with key-value pairs
                and `SAMPLES` are a list of values for each sample. Each value
                is also a list of values for each FORMAT.
            If a record matches multiple fixes, the first one is applied.
        helpers: raw code the provide some helpers for the fixes
            The code will automatically dedented if given as a string. A list
            of strings is also supported and will be joined with newlines.
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.python
    envs = {"fixes": [], "helpers": ""}
    script = "file://../scripts/vcf/VcfFix.py"


class VcfAnno(Proc):
    """Annotate a VCF file using vcfanno

    https://github.com/brentp/vcfanno

    Input:
        infile: The input VCF file
        conffile: The configuration file for vcfanno or configuration dict
            itself

    Output:
        outfile: The output VCF file

    Envs:
        vcfanno: Path to vcfanno
        ncores: Number of cores to use
        conffile: configuration file for vcfanno or configuration dict itself
            This is ignored when `conffile` is given as input
        args: Additional arguments to pass to vcfanno

    Requires:
        - name: vcfanno
          check: |
            {{proc.envs.vcfanno}} --help
    """

    input = "infile:file, conffile"
    output = "outfile:file:{{in.infile | stem0}}.{{envs.tool}}.vcf"
    lang = config.lang.python
    envs = {
        "vcfanno": config.exe.vcfanno,
        "ncores": config.misc.ncores,
        "conffile": {},
        "args": {"permissive-overlap": True},
    }
    script = "file://../scripts/vcf/VcfAnno.py"


class TruvariBench(Proc):
    """Run `truvari bench` to compare a VCF with CNV calls and
    base CNV standards

    Requires truvari v4+

    See https://github.com/ACEnglish/truvari/wiki/bench

    Input:
        compvcf: The VCF file with CNV calls to compare
        basevcf: The VCF file with standard CNVs

    Output:
        outdir: The output directory

    Envs:
        truvari: Path to truvari
        `<other>`: Ohter `truvari bench` arguments

    Requires:
        truvari:
            - check: {{proc.envs.truvari}} version
    """
    input = "compvcf:file, basevcf:file"
    output = "outdir:dir:{{in.compvcf | stem0 | append: '.truvari_bench'}}"
    envs = {
        "truvari": config.exe.truvari,
        "ref": config.ref.reffa,
        "refdist": 500,
        "pctseq": 0.7,
        "pctsize": 0.7,
        "pctovl": 0.0,
        "typeignore": False,
        "multimatch": False,
    }
    lang = config.lang.bash
    script = "file://../scripts/vcf/TruvariBench.sh"


class TruvariBenchSummary(Proc):
    """Summarise the statistics from `TruvariBench` for multiple jobs (VCFs)

    Input:
        indirs: The input directories, which should be the output directories
            of `TruvariBench`

    Output:
        outdir: The output directory, including the summary table and plots

    Envs:
        plots: The stats to plot with barplots.
            Candidates are `TP-base`, `TP-call`, `FP`, `FN`, `precision`,
            `recall`, `f1`, `base cnt`, `call cnt`, `TP-call_TP-gt`,
            `TP-call_FP-gt`, `TP-base_TP-gt`, `TP-base_FP-gt`, and
            `gt_concordance`
            See https://github.com/ACEnglish/truvari/wiki/bench
        devpars: The parameters to use for the plots.

    Requires:
        r-ggprism:
            - check: {{proc.lang}} -e "library(ggprism)"
        r-rjson:
            - check: {{proc.lang}} -e "library(rjson)"
        r-dplyr:
            - check: {{proc.lang}} -e "library(dplyr)"
        r-ggplot2:
            - check: {{proc.lang}} -e "library(ggplot2)"
    """
    input = "indirs:files"
    input_data = lambda ch: [list(ch.iloc[:, 0])]
    output = "outdir:dir:truvari_bench.summary"
    lang = config.lang.rscript
    envs = {
        "plots": ["comp cnt", "base cnt", "precision", "recall", "f1"],
        "devpars": None,
    }
    script = "file://../scripts/vcf/TruvariBenchSummary.R"
    plugin_opts = {"report": "file://../reports/vcf/TruvariBenchSummary.svelte"}


class TruvariConsistency(Proc):
    """Run `truvari consistency` to check consistency of CNV calls

    See https://github.com/ACEnglish/truvari/wiki/consistency

    Requires truvari v4+

    Input:
        vcfs: The vcf files with CNV calls

    Output:
        outfile: The output file with the report

    Envs:
        truvari: Path to truvari
        heatmap: Whether to generate a heatmap of the consistency
            Set to False to disable
            annofile: The annotation file for the heatmap, multiple columns
            but the first column must be the sample name. Note that the stem of
            the vcf file name from consistency file will be used. These
            annotations will be added as row annotations.
            Other options see also `biopipen.ns.plot.Heatmap`.
    """
    input = "vcfs:files"
    output = (
        "outdir:dir:"
        "{{in.vcfs | first | stem0 | append: '.etc.truvari_consistency'}}"
    )
    lang = config.lang.rscript
    envs = {"truvari": config.exe.truvari, "heatmap": {}}
    script = "file://../scripts/vcf/TruvariConsistency.R"
    plugin_opts = {"report": "file://../reports/vcf/TruvariConsistency.svelte"}


class BcftoolsAnnotate(Proc):
    """Add or remove annotations from VCF files

    See also: <https://samtools.github.io/bcftools/bcftools.html#annotate>

    Input:
        infile: The input VCF file
        annfile: The annotation file.
            Currently only VCF files are supported.

    Output:
        outfile: The VCF file with annotations added or removed.

    Envs:
        bcftools: Path to bcftools
        tabix: Path to tabix, used to index infile and annfile
        annfile: The annotation file. If `in.annfile` is provided,
            this is ignored
        ncores (type=int): Number of cores (`--threads`) to use
        columns (auto): Comma-separated or list of columns or tags to carry over from
            the annotation file. Overrides `-c, --columns`
        remove (auto): Remove the specified columns from the input file
        header (list): Headers to be added
        gz (flag): Whether to gzip the output file
        index (flag): Whether to index the output file (tbi) (`envs.gz` forced to True)
        <more>: Other arguments for `bcftools annotate`
            See also <https://samtools.github.io/bcftools/bcftools.html#annotate>
            Note that the underscore `_` will be replaced with dash `-` in the
            argument name.
    """
    input = "infile:file, annfile:file"
    output = (
        "outfile:file:{{in.infile | stem: 'gz'}}.vcf"
        "{{'.gz' if envs.index or envs.gz else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "annfile": None,
        "columns": [],
        "remove": [],
        "header": [],
        "gz": True,
        "index": True,
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/vcf/BcftoolsAnnotate.py"


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
    script = "file://../scripts/vcf/BcftoolsFilter.py"


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
    script = "file://../scripts/vcf/BcftoolsSort.py"


class BcftoolsMerge(Proc):
    """Merge multiple VCF files using `bcftools merge`.

    Input:
        infiles: The input VCF files

    Output:
        outfile: The merged VCF file.

    Envs:
        bcftools: Path to bcftools
        tabix: Path to tabix, used to index infile/outfile
        ncores (type=int): Number of cores (`--threads`) to use
        gz (flag): Whether to gzip the output file
        index (flag): Whether to index the output file (tbi) (`envs.gz` forced to True)
        <more>: Other arguments for `bcftools merge`.
            See also <https://samtools.github.io/bcftools/bcftools.html#merge>
    """
    input = "infiles:files"
    output = (
        "outfile:file:{{in.infiles | first | stem | append: '_etc_merged'}}.vcf"
        "{{'.gz' if envs.index or envs.gz else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
        "gz": True,
        "index": True,
    }
    script = "file://../scripts/vcf/BcftoolsMerge.py"


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
    script = "file://../scripts/vcf/BcftoolsView.py"
