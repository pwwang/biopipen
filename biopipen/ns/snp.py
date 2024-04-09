"""Plink processes"""

from ..core.proc import Proc
from ..core.config import config


class PlinkSimulation(Proc):
    """Simulate SNPs using PLINK v1.9

    See also <https://www.cog-genomics.org/plink/1.9/input#simulate> and
    <https://pwwang.github.io/biopipen/api/biopipen.ns.snp/#biopipen.ns.snp.PlinkSimulation>

    Input:
        configfile: Configuration file containing the parameters for the simulation.
            The configuration file (in toml, yaml or json format) should contain a
            dictionary of parameters.  The parameters are listed in `envs` except
            `ncores`, which is used for parallelization. You can set parameters
            in `envs` and override them in the configuration file.

    Output:
        outdir: Output directory containing the simulated data
            `plink_sim.bed`, `plink_sim.bim`, and `plink_sim.fam` will be generated.
        gtmat: Genotype matrix file containing the simulated data with rows representing
            SNPs and columns representing samples.

    Envs:
        nsnps (type=int): Number of SNPs to simulate
        ncases (type=int): Number of cases to simulate
        nctrls (type=int): Number of controls to simulate
        plink: Path to PLINK v1.9
        seed (type=int): Random seed. If not set, seed will not be set.
        label: Prefix label for the SNPs.
        prevalence  (type=float): Disease prevalence.
        minfreq (type=float): Minimum allele frequency.
        maxfreq (type=float): Maximum allele frequency.
        hetodds (type=float): Odds ratio for heterozygous genotypes.
        homodds (type=float): Odds ratio for homozygous genotypes.
        missing (type=float): Proportion of missing genotypes.
        args (ns): Additional arguments to pass to PLINK.
            - <more>: see <https://www.cog-genomics.org/plink/1.9/input#simulate>.
        transpose_gtmat (flag): If set, the genotype matrix (`out.gtmat`) will
            be transposed.
        sample_prefix: Use this prefix for the sample names. If not set, the sample
            names will be `per0_per0`, `per1_per1`, `per2_per2`, etc. If set, the
            sample names will be `prefix0`, `prefix1`, `prefix2`, etc.
            This only affects the sample names in the genotype matrix file
            (`out.gtmat`).
    """
    input = "configfile:file"
    output = [
        "outdir:dir:{{in.configfile | stem}}.plink_sim",
        "gtmat:file:{{in.configfile | stem}}.plink_sim/"
        "{{in.configfile | stem}}-gtmat.txt",
    ]
    lang = config.lang.python
    envs = {
        "nsnps": None,
        "ncases": None,
        "nctrls": None,
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
        "transpose_gtmat": False,
        "sample_prefix": None,
    }
    script = "file://../scripts/snp/PlinkSimulation.py"


class MatrixEQTL(Proc):
    """Run Matrix eQTL

    See also <https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/>

    Input:
        geno: Genotype matrix file with rows representing SNPs and columns
            representing samples.
        expr: Expression matrix file with rows representing genes and columns
            representing samples.
        cov: Covariate matrix file with rows representing covariates and columns
            representing samples.

    Output:
        alleqtls: Matrix eQTL output file
        cisqtls: The cis-eQTL file if `snppos` and `genepos` are provided.
            Otherwise it'll be empty.

    Envs:
        model (choice): The model to use.
            - `linear`: Linear model
            - `modelLINEAR`: Same as `linear`
            - `anova`: ANOVA model
            - `modelANOVA`: Same as `anova`
        pval (type=float): P-value threshold for eQTLs
        transp (type=float): P-value threshold for trans-eQTLs.
            If cis-eQTLs are not enabled (`snppos` and `genepos` are not set),
            this defaults to 1e-5.
            If cis-eQTLs are enabled, this defaults to `None`, which will disable
            trans-eQTL analysis.
        fdr (flag): Do FDR calculation or not (save memory if not).
        snppos: The path of the SNP position file.
            It could be a BED, GFF, VCF or a tab-delimited file with
            `snp`, `chr`, `pos` as the first 3 columns.
        genepos: The path of the gene position file.
            It could be a BED or GFF file.
        dist (type=int): Distance threshold for cis-eQTLs.
        transpose_geno (flag): If set, the genotype matrix (`in.geno`)
            will be transposed.
        transpose_expr (flag): If set, the expression matrix (`in.expr`)
            will be transposed.
        transpose_cov (flag): If set, the covariate matrix (`in.cov`)
            will be transposed.
    """
    input = "geno:file, expr:file, cov:file"
    output = [
        "alleqtls:file:{{in.geno | stem}}.alleqtls.txt",
        "cisqtls:file:{{in.geno | stem}}.cisqtls.txt",
    ]
    lang = config.lang.rscript
    envs = {
        "model": "linear",
        "pval": 1e-3,
        "transp": None,
        "fdr": False,
        "snppos": None,
        "genepos": config.ref.refgene,
        "dist": 250000,
        "transpose_geno": False,
        "transpose_expr": False,
        "transpose_cov": False,
    }
    script = "file://../scripts/snp/MatrixEQTL.R"
