"""Plink processes"""

from ..core.proc import Proc
from ..core.config import config


class PlinkSimulation(Proc):
    """Simulate SNPs using PLINK v1.9

    See also <https://www.cog-genomics.org/plink/1.9/input#simulate>.

    Input:
        nsnps: Number of SNPs to simulate
        ncases: Number of cases to simulate
        nctrls: Number of controls to simulate

    Output:
        outdir: Output directory containing the simulated data
            `plink_sim.bed`, `plink_sim.bim`, and `plink_sim.fam` will be generated.

    Envs:
        plink: Path to PLINK v1.9
        seed (type=int): Random seed.
            If not set, seed will not be set.
        label: Prefix label for the SNPs.
        prevalence  (type=float): Disease prevalence.
        minfreq (type=float): Minimum allele frequency.
        maxfreq (type=float): Maximum allele frequency.
        hetodds (type=float): Odds ratio for heterozygous genotypes.
        homodds (type=float): Odds ratio for homozygous genotypes.
        missing (type=float): Proportion of missing genotypes.
        args (ns): Additional arguments to pass to PLINK.
            - <more>: see <https://www.cog-genomics.org/plink/1.9/input#simulate>.
    """
    input = "nsnps:var, ncases:var, nctrls:var"
    output = (
        "outdir:dir:{{in.nsnps | int}}_"
        "{{in.ncases | int}}xcases_{{in.nctrls | int}}xctrls.plink_sim"
    )
    lang = config.lang.python
    envs = {
        "plink": config.exe.plink,
        "seed": None,
        "label": "SNP",
        "prevalence": 0.01,
        "minfreq": 0.0,
        "maxfreq": 1.0,
        "hetodds": 1.0,
        "homodds": 1.0,
        "missing": 0.0,
        "args": {},
    }
    script = "file://../scripts/snp/PlinkSimulation.py"
