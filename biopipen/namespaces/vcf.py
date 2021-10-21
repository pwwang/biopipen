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
