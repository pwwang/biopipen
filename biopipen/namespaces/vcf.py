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
    envs = {
        "tabix": config.exe.tabix,
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
