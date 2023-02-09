"""Processes for TCGA MAF files."""
from ..core.proc import Proc
from ..core.config import config


class Maf2Vcf(Proc):
    """Converts a MAF file to a VCF file.

    This is a wrapper around the `maf2vcf` script from the `maf2vcf` package.

    Input:
        infile: The input MAF file

    Output:
        outfile: Output multi-sample VCF containing all TN-pairs
        outdir: Path to output directory where VCFs will be stored,
            one per TN-pair

    Envs:
        perl: Path to perl to run `maf2vcf.pl`
        samtools: Path to samtools to be used in `maf2vcf.pl`
        args: Other arguments to pass to the script
    """
    input = "infile:file"
    output = [
        'outfile:file:{{in.infile | stem}}.vcfs/'
        '{{in.infile | stem}}.multisample.vcf',
        'outdir:dir:{{in.infile | stem}}.vcfs'
    ]
    lang = config.lang.python
    envs = {
        "perl": config.lang.perl,
        "samtools": config.exe.samtools,
        "ref": config.ref.reffa,
        "args": {"per-tn-vcfs": True},
    }
    script = "file://../scripts/tcgamaf/Maf2Vcf.py"


class MafAddChr(Proc):
    """Adds the `chr` prefix to chromosome names in a MAF file if not present.

    Input:
        infile: The input MAF file

    Output:
        outfile: The output MAF file
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.maf"
    lang = config.lang.python
    script = "file://../scripts/tcgamaf/MafAddChr.py"
