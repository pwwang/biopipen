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

    Envs:
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
            Any converts return `None` will skip the record
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
        "genome": config.ref.genome,
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


class BedConsensus(Proc):
    """Find consensus regions from multiple BED files.

    Unlike `bedtools merge/cluster`, it does not find the union regions nor
    intersect regions. Instead, it finds the consensus regions using the
    distributions of the scores of the bins

                                         bedtools cluster
    Bedfile A            |----------|    1
    Bedfile B          |--------|        1
    Bedfile C              |------|      1
    BedConsensus         |--------|
    bedtools intesect      |----|
    bedtools merge     |------------|
    Distribution       |1|2|3333|2|1|    (later normalized into 0~1)

    If column #5 is provided, it can be used as weights to determine the
    ends for the consensus regions. The weight score should be the sum of
    all base pairs in the region.

    Input:
        bedfiles: Input BED files

    Output:
        outbed: The output BED file

    Envs:
        bedtools: The path to bedtools
        binsize: The binsize to calculate the weights
            The smaller the more accurate to determine the ends of
            consensus regions.
        ignore_scores: A list of indices of the BED files to ignore their
            scores (column #5)
            If scores are ignored or not provided, use `1.0`.
        cutoff: The cutoff to determine the ends of consensus regions
            The cutoff weights of bins are used to determine the ends
        distance: When the distance between two bins is smaller than this value,
            they are merged into one bin using `bedtools merge -d`. `0` means
            no merging.
        chrsize: The chromosome size file, used to make windows
        ncores: Number of cores to use to calculate the weights for
            each bed file
    """
    input = "bedfiles:files"
    output = (
        "outbed:file:{{in.bedfiles | first | stem | append: '_consensus'}}.bed"
    )
    lang = config.lang.python
    envs = {
        "bedtools": config.exe.bedtools,
        "binsize": 1000,
        "ncores": config.misc.ncores,
        "ignore_scores": [],
        "cutoff": 0.5,
        "distance": 1,
        "chrsize": config.ref.chrsize,
    }
    script = "file://../scripts/bed/BedConsensus.py"
