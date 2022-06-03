"""Tools to handle BED files"""
from ..core.proc import Proc
from ..core.config import config

class BedLiftOver(Proc):
    """Liftover a BED file using liftOver

    Input:
        inbed: The input BED file

    Output:
        outbed: The output BED file

    Envs:
        liftover: The path to liftOver
        chain: The map chain file for liftover

    Requires:
        - name: liftOver
          check: |
            {{proc.envs.liftover}} 2>&1 | grep "usage"
    """
    input = "inbed:file"
    output = "outbed:file:{{in.inbed | basename}}"
    envs = {
        "liftover": config.exe.liftover,
        "chain": config.path.liftover_chain,
    }
    lang = config.lang.bash
    script = "file://../scripts/bed/BedLiftOver.sh"



class Bed2Vcf(Proc):
    """Convert a BED file to a valid VCF file with minimal information

    Input:
        inbed: The input BED file

    Output:
        outvcf: The output VCF file

    Args:
        sample: The sample name to be used in the VCF file
            You can use a lambda function (in string) to generate
            the sample name from the stem of input file
        ref: The reference fasta file, used to grab the reference allele.
            To add contigs in header, the `fai` file is also required at
            `<ref>.fai`
        genome: The genome assembly, added as `source` in header
        base: 0 or 1, whether the coordinates in BED file are 0- or 1-based
        headers: The header lines to be added to the VCF file
        infos: The INFO dicts to be added to the VCF file
        formats: The FORMAT dicts to be added to the VCF file
            The keys 'ID', 'Description', 'Type', and 'Number' are required.
        converters: A dict of converters to be used for each INFO or FORMAT
            The key is the ID of an INFO or FORMAT, and the value is
        nonexisting_contigs: Whether to `keep` or `drop` the non-existing
            contigs in `ref`.
        helpers: Raw code to be executed to provide some helper functions
            since only lambda functions are supported in converters
        index: Sort and index output file

    Requires:
        - name: cyvcf2
          check: |
            {{proc.lang}} -c "import cyvcf2"
        - name: pysam
          check: |
            {{proc.lang}} -c "import pysam"
        - name: bcftools
          if: {{proc.envs.index}}
          check: |
            {{proc.envs.bcftools}} --version
    """
    input = "inbed:file"
    output = (
        "outvcf:file:{{in.inbed | stem}}.vcf{{'.gz' if envs.index else ''}}"
    )
    lang = config.lang.python
    envs = {
        "bcftools": config.exe.bcftools,
        "sample": "lambda stem: stem",
        "ref": config.ref.reffa,
        "genome": "",
        "nonexisting_contigs": "drop",
        "base": 0,
        "index": True,
        "headers": [],
        "infos": [],
        "formats": [],
        "converters": {},
        "helpers": "",
    }
    script = "file://../scripts/bed/Bed2Vcf.py"