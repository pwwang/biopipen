"""Tools to handle VCF files"""
from ..core.proc import Proc
from ..core.params import params

class VcfLiftOver(Proc):
    """Liftover a VCF file using GATK

    Input:
        invcf: The input VCF file

    Output:
        outvcf: The output VCF file

    Envs:
        chain: The map chain file for liftOver
    """
    input = "invcf:file"
    output = "outvcf:file:{{in.invcf | stem}}"
    envs = {
        "liftover": params.exe.liftover,
        "chain": params.path.liftover_chain,
    }
    lang = params.lang.bash
    script = "file://../scripts/vcf/VcfLiftOver.sh"
