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
    """

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
    script = "file://../scripts/vcf/VcfDownSample.sh"


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

    Requires:
        - name: biopipen
          check: |
            {{proc.lang}} -c "import biopipen"

    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.python
    envs = {"fixes": [], "helpers": ""}
    script = "file://../scripts/vcf/VcfFix.py"


class TruvariBench(Proc):
    """Run `truvari bench` to compare a VCF with CNV calls and
    base CNV standards

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
        - name: truvari
          check: |
            {{proc.envs.truvari}} version
    """
    input = "compvcf:file, basevcf:file"
    output = "outdir:dir:{{in.compvcf | stem0 | append: '.truvari_bench'}}"
    envs = {
        "truvari": config.exe.truvari,
        "ref": config.ref.reffa,
        "refdist": 500,
        "pctsim": 0.7,
        "pctsize": 0.7,
        "pctovl": 0.0,
        "typeignore": False,
    }
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
        - name: r-ggprism
          check: |
            {{proc.lang}} -e "library(ggprism)"
        - name: r-rjson
          check: |
            {{proc.lang}} -e "library(rjson)"
        - name: r-dplyr
          check: |
            {{proc.lang}} -e "library(dplyr)"
        - name: r-ggplot2
          check: |
            {{proc.lang}} -e "library(ggplot2)"

    """
    input = "indirs:files"
    input_data = lambda ch: [list(ch.iloc[:, 0])]
    output = "outdir:dir:truvari_bench.summary"
    lang = config.lang.rscript
    envs = {
        "plots": ["call cnt", "base cnt", "precision", "recall", "f1"],
        "devpars": None,
    }
    script = "file://../scripts/vcf/TruvariBenchSummary.R"
    plugin_opts = {
        "report": "file://../reports/vcf/TruvariBenchSummary.svelte"
    }


class TruvariConsistency(Proc):
    """Run `truvari consistency` to check consistency of CNV calls

    See https://github.com/ACEnglish/truvari/wiki/consistency

    Input:
        vcfs: The vcf files with CNV calls

    Output:
        outfile: The output file with the report

    Envs:
        truvari: Path to truvari
    """
    input = "vcfs:files"
    output = (
        "outfile:file:"
        "{{in.vcfs | first | stem0 | append: '.etc.truvari_consistency.txt'}}"
    )
    envs = {"truvari": config.exe.truvari}
    script = "file://../scripts/vcf/TruvariConsistency.sh"
