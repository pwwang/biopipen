"""Plink processes"""

from ..core.proc import Proc
from ..core.config import config


class PlinkSimulation(Proc):
    """Simulate SNPs using PLINK v2

    See also <https://www.cog-genomics.org/plink/2.0/input#simulate> and
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
        plink: Path to PLINK v2
        seed (type=int): Random seed. If not set, seed will not be set.
        label: Prefix label for the SNPs.
        prevalence  (type=float): Disease prevalence.
        minfreq (type=float): Minimum allele frequency.
        maxfreq (type=float): Maximum allele frequency.
        hetodds (type=float): Odds ratio for heterozygous genotypes.
        homodds (type=float): Odds ratio for homozygous genotypes.
        missing (type=float): Proportion of missing genotypes.
        args (ns): Additional arguments to pass to PLINK.
            - <more>: see <https://www.cog-genomics.org/plink/2.0/input#simulate>.
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
            - linear: Linear model
            - modelLINEAR: Same as `linear`
            - anova: ANOVA model
            - modelANOVA: Same as `anova`
        pval (type=float): P-value threshold for eQTLs
        match_samples (flag): Match samples in the genotype and expression matrices.
            If True, an error will be raised if samples from `in.geno`, `in.expr`,
            and `in.cov` (if provided) are not the same.
            If False, common samples will be used to subset the matrices.
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
        "match_samples": False,
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


class PlinkFromVcf(Proc):
    """Convert VCF to PLINK format.

    The PLINK format consists of 3 files: `.bed`, `.bim`, and `.fam`.

    Requires PLINK v2

    TODO:
        Handle sex when sex chromosomes are included.

    Input:
        invcf: VCF file

    Output:
        outdir: Output directory containing the PLINK files

    Envs:
        plink: Path to PLINK v2
        tabix: Path to tabix
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        vcf_half_call (choice): The current VCF standard does not specify
            how '0/.' and similar GT values should be interpreted.
            - error: error out and reports the line number of the anomaly
            - e: alias for `error`
            - haploid: treat half-calls as haploid/homozygous
            - h: alias for `haploid`
            - missing: treat half-calls as missing
            - m: alias for `missing`
            - reference: treat the missing part as reference
            - r: alias for `reference`
        double_id (flag): set both FIDs and IIDs to the VCF/BCF sample ID.
        vcf_filter (auto): skip variants which failed one or more filters tracked
            by the FILTER field.
            If True, only FILTER with `PASS` or `.` will be kept.
            Multiple filters can be specified by separating them with space or
            as a list.
        vcf_idspace_to: convert all spaces in sample IDs to this character.
        set_missing_var_ids: update variant IDs using a template string,
            with a '@' where the chromosome code should go, and a '#' where the
            base-pair position belongs. You can also specify `\\$r` and `\\$a` for
            the reference and alternate alleles, respectively.
            See <https://www.cog-genomics.org/plink/2.0/data#set_all_var_ids>
        max_alleles (type=int): Maximum number of alleles per variant.
        <more>: see <https://www.cog-genomics.org/plink/2.0/> for more options.
            Note that `_` will be replaced by `-` in the argument names.
    """  # noqa: E501
    input = "invcf:file"
    output = "outdir:dir:{{in.invcf.stem | regex_replace: '\\.gz$', ''}}"
    lang = config.lang.python
    envs = {
        "plink": config.exe.plink2,
        "tabix": config.exe.tabix,
        "ncores": config.misc.ncores,
        "vcf_half_call": "missing",
        "double_id": True,
        "vcf_filter": True,
        "vcf_idspace_to": "_",
        "set_missing_var_ids": "@_#",
        "max_alleles": 2,
    }
    script = "file://../scripts/snp/PlinkFromVcf.py"


class Plink2GTMat(Proc):
    """Convert PLINK files to genotype matrix.

    Requires PLINK v2. The .raw/.traw file is generated by plink and then transformed
    to a genotype matrix file.
    See <https://www.cog-genomics.org/plink/2.0/formats#raw> and
    <https://www.cog-genomics.org/plink/2.0/formats#traw> for more information.

    The allelic dosage is used as the values of genotype matrix.
    "--keep-allele-order" is used to keep the allele order consistent with the
    reference allele first. This way, the genotype of homozygous reference alleles
    will be encoded as 2, heterozygous as 1, and homozygous alternate alleles as 0.
    This is the PLINK dosage encoding. If you want to use this encoding, you can
    set `envs.gtcoding` to `plink`. Otherwise, the default encoding is `vcf`, which
    will encode the genotype as 0, 1, and 2 for homozygous reference, heterozygous,
    and homozygous alternate alleles, respectively.

    Note that `envs.gtcoding = "vcf"` only works for biallelic variants for now.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outfile: Genotype matrix file with rows representing SNPs and columns
            representing samples if `envs.transpose` is `False`.

    Envs:
        plink: Path to PLINK v2.0
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        transpose (flag): If set, the genotype matrix (`out.outfile`) is transposed.
        samid: what to use as sample ID.
            Placeholders include `{fid}` and `{iid}` for family and individual IDs,
            respectively.
        varid: what to use as variant ID.
            Placeholders include `{chr}`, `{pos}`, `{rs}`, `{ref}`, and `{alt}` for
            chromosome, position, rsID, reference allele, and alternate allele,
            respectively.
        trans_chr: A dictionary to translate chromosome numbers to chromosome names.
        missing_id: what to use as the rs if missing.
        gtcoding (choice): The genotype coding to use.
            - vcf: 0/1/2 for homozygous reference, heterozygous, and homozygous
                alternate alleles, respectively.
            - plink: 2/1/0 for homozygous reference, heterozygous, and homozygous
                alternate alleles, respectively.
    """
    input = "indir:dir"
    output = "outfile:file:{{in.indir | stem}}-gtmat.txt"
    lang = config.lang.python
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "transpose": False,
        "samid": "{fid}_{iid}",
        "varid": "{chr}_{pos}_{varid}_{ref}_{alt}",
        "trans_chr": {"23": "X", "24": "Y", "25": "XY", "26": "M"},
        "missing_id": "NA",
        "gtcoding": "vcf",
    }
    script = "file://../scripts/snp/Plink2GTMat.py"


class PlinkIBD(Proc):
    """Run PLINK IBD analysis (identity by descent)

    See also <https://www.cog-genomics.org/plink/1.9/ibd>
    This has to run with PLINK v1.9. Plink v2 does not support IBD analysis yet.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outdir: Output file containing the IBD results.
            Including [`.genome`](https://www.cog-genomics.org/plink/2.0/formats#genome)
            file for the original IBD report from PLINK, and `.ibd.png` for the
            heatmap of `PI_HAT` values.

    Envs:
        plink: Path to PLINK v1.9
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        highld: High LD regions to be excluded from the analysis.
            If not set, no regions will be excluded.
        samid: what to use as sample ID.
            Placeholders include `{fid}` and `{iid}` for family and individual IDs,
            respectively
        indep (type=auto): LD pruning parameters. Either a list of numerics or a string
            concatenated by `,` to specify
            1) consider a window of N SNPs (e.g. 50),
            2) calculate LD between each pair of SNPs in the window (e.g. 5),
            3) remove one of a pair of SNPs if the LD is greater than X (e.g. 0.2).
        pihat (type=float): PI_HAT threshold for IBD analysis.
            See also <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007749/>
        plot (flag): If set, plot the heatmap of `PI_HAT` values.
        anno: The annotation file for the samples, used to plot on the heatmap.
            Names must match the ones that are transformed by `args.samid`.
        seed (type=int): Random seed for the analysis.
        devpars (ns): The device parameters for the plot.
            - width (type=int): Width of the plot
            - height (type=int): Height of the plot
            - res (type=int): Resolution of the plot
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}.ibd"
    lang = config.lang.rscript
    envs = {
        "plink": config.exe.plink,
        "ncores": config.misc.ncores,
        "highld": None,
        "samid": "{fid}_{iid}",
        "indep": [50, 5, 0.2],
        "pihat": 0.1875,
        "plot": True,
        "anno": None,
        "seed": 8525,
        "devpars": {"width": 1000, "height": 1000, "res": 100},
    }
    script = "file://../scripts/snp/PlinkIBD.R"
    plugin_opts = {"report": "file://../reports/snp/PlinkIBD.svelte"}


class PlinkHWE(Proc):
    """Hardy-Weinberg Equilibrium report and filtering

    See also <https://www.cog-genomics.org/plink/2.0/basic_stats#hardy>

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outdir: Output file containing the HWE results.
            Including [`.hwe`](https://www.cog-genomics.org/plink/2.0/formats#hwe)
            file for the original HWE report from PLINK and
            `.hardy.fail` for the variants that failed the HWE test.
            It also includes binary files `.bed`, `.bim`, and `.fam`

    Envs:
        plink: Path to PLINK v2
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        cutoff (type=float): P-value cutoff for HWE test
        plot (flag): If set, plot the distribution of HWE p-values.
        devpars (ns): The device parameters for the plot.
            - width (type=int): Width of the plot
            - height (type=int): Height of the plot
            - res (type=int): Resolution of the plot
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}.hwe"
    lang = config.lang.rscript
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "cutoff": 1e-5,
        "plot": True,
        "devpars": {"width": 1000, "height": 800, "res": 100},
    }
    script = "file://../scripts/snp/PlinkHWE.R"
    plugin_opts = {"report": "file://../reports/snp/PlinkHWE.svelte"}


class PlinkHet(Proc):
    """Calculation of sample heterozygosity.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outdir: Output file containing the heterozygosity results.
            Including [`.het`](https://www.cog-genomics.org/plink/2.0/formats#het)
            file for the original heterozygosity report from PLINK and
            `.het.fail` for the samples that failed the heterozygosity test.
            It also includes binary files `.bed`, `.bim`, and `.fam`

    Envs:
        plink: Path to PLINK v2, at least v2.00a5.10
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        cutoff (type=float): Heterozygosity cutoff, samples with heterozygosity
            beyond `mean - cutoff * sd` or `mean + cutoff * sd` will be considered
            as outliers.
        plot (flag): If set, plot the distribution of heterozygosity values.
        devpars (ns): The device parameters for the plot.
            - width (type=int): Width of the plot
            - height (type=int): Height of the plot
            - res (type=int): Resolution of the plot
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}.het"
    lang = config.lang.rscript
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "cutoff": 3.0,
        "plot": True,
        "devpars": {"width": 1000, "height": 800, "res": 100},
    }
    script = "file://../scripts/snp/PlinkHet.R"
    plugin_opts = {"report": "file://../reports/snp/PlinkHet.svelte"}


class PlinkCallRate(Proc):
    """Calculation of call rate for the samples and variants.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outdir: Output file containing the call rate results.
            Including [`.imiss`](https://www.cog-genomics.org/plink/2.0/formats#imiss)
            file for missing calls for samples,
            [`.lmiss`](https://www.cog-genomics.org/plink/2.0/formats#lmiss) for
            missing calls for variants, `.samplecr.fail` for the samples fail
            sample call rate cutoff (`args.samplecr`), and `.varcr.fail` for the SNPs
            fail snp call rate cutoff (`args.varcr`).
            It also includes binary files `.bed`, `.bim`, and `.fam`.

    Envs:
        plink: Path to PLINK v2
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        samplecr (type=float): Sample call rate cutoff
        varcr (type=float): Variant call rate cutoff
        max_iter (type=int): Maximum number of iterations to run the call rate
            calculation.
            Since the sample and variant call rates are affected by each other,
            it may be necessary to iterate the calculation to get the stable results.
        plot (flag): If set, plot the distribution of call rates.
        devpars (ns): The device parameters for the plot.
            - width (type=int): Width of the plot
            - height (type=int): Height of the plot
            - res (type=int): Resolution of the plot
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}.callrate"
    lang = config.lang.rscript
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "samplecr": 0.95,
        "varcr": 0.95,
        "max_iter": 3,
        "plot": True,
        "devpars": {"width": 1000, "height": 800, "res": 100},
    }
    script = "file://../scripts/snp/PlinkCallRate.R"
    plugin_opts = {"report": "file://../reports/snp/PlinkCallRate.svelte"}


class PlinkFilter(Proc):
    """Filter samples and variants for PLINK files.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files
        samples_file: File containing the sample IDs.
        variants_file: File containing the variant IDs or regions.

    Output:
        outdir: Output directory containing the filtered PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Envs:
        plink: Path to PLINK v2
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        samples (auto): Sample IDs.
            If both FID and IID should be provided and separatedby `/`. Otherwise,
            assuming the same FID and IID.
            A list of sample IDs or string concatenated by `,`.
            If either `in.samples_file` or `envs.samples_file` is set,
            this will be ignored.
        variants (auto): Variant IDs.
            A list of variant IDs or string concatenated by `,`.
            If either `in.variants_file` or `envs.variants_file` is set,
            this will be ignored.
        samples_file: File containing the sample IDs.
            If `in.samples_file` is set, this will be ignored.
        variants_file: File containing the variant IDs.
            If `in.variants_file` is set, this will be ignored.
        keep (flag): Use `samples`/`variants`/`samples_file`/`variants_file` to
            only keep the specified samples/variants, instead of removing them.
        vfile_type (choice): The type of the variants file.
            - id: Variant IDs
            - bed0: 0-based BED file
            - bed1: 1-based BED file
        chr: Chromosome to keep.
            For example, `1-4 22 XY` will keep chromosomes 1 to 4, 22, and XY.
        not_chr: Chromosome to remove.
            For example, `1-4 22 XY` will remove chromosomes 1 to 4, 22, and XY.
        autosome (flag): Excludes all unplaced and non-autosomal variants
        autosome_xy (flag): Does `autosome` but does not exclude the pseudo-autosomal
            region of X.
        snps_only (auto): Excludes all variants with one or more multi-character
            allele codes. With 'just-acgt', variants with single-character allele codes
            outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', <missing code>}
            are also excluded.
    """
    input = [
        "indir:dir",
        "samples_file:file",
        "variants_file:file",
    ]
    output = "outdir:dir:{{in.indir | stem}}.filtered"
    lang = config.lang.python
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "samples": None,
        "variants": None,
        "samples_file": None,
        "variants_file": None,
        "keep": False,
        "vfile_type": "id",
        "chr": None,
        "not_chr": None,
        "autosome": False,
        "autosome_xy": False,
        "snps_only": False,
    }
    script = "file://../scripts/snp/PlinkFilter.py"


class PlinkFreq(Proc):
    """Calculate allele frequencies for the variants.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Output:
        outdir: Output file containing the allele frequency results.
            By default, it includes
            [`.afreq`](https://www.cog-genomics.org/plink/2.0/formats#afreq)
            file for the allele frequency report from PLINK.
            Modifiers can be added to change this behavior.
            See `envs.modifier` for more information.
            When `envs.filter != no`, it also includes binary files `.bed`, `.bim`,
            and `.fam` after filtering with `envs.cutoff`.

    Envs:
        plink: Path to PLINK v2
        ncores (type=int): Number of cores/threads to use, will pass to plink
            `--threads` option
        modifier (choice): The modifier of `--freq` to control the output behavior.
            - none: No modifier, only the `.afreq` file will be generated.
                `MAF` (minor allele frequency) will be added in addition to the
                `REF_FREQ` and `ALT1_FREQ` columns. Check `.afreqx` for the added
                columns.
            - counts: write allele count report to `.acount`.
                See <https://www.cog-genomics.org/plink/2.0/formats#afreq>.
                `ALT1`, `ALT1_CT`, and `REF_CT` are added. Check `.acountx` for
                the added columns.
            - x: write genotype count report to `.gcount`
                Like `--freqx` in v1.9, `--geno-counts` will be run to generate
                the genotype counts.
                `ALT1`, `HET_REF_ALT1_CT`, and `HOM_ALT1_CT` are added. Check
                `.gcountx` for the added columns.
        gz (flag): If set, compress the output files.
        cutoff (auto): Cutoffs to mark or filter the variants.
            If a float is given, default column will be used based on the modifier.
            For `modifier="none"`, it defaults to `MAF`.
            For `modifier="counts"`, it defaults to `ALT1_CT`.
            For `modifier="x"`, it defaults to `HOM_ALT1_CT`.
            Or this could be a dictionary to specify the column names and cutoffs.
            For example, `{"MAF": 0.05}`.
        filter (auto): The direction of filtering variants based on `cutoff`.
            If a single value is given, it will apply to all columns provided in
            `cutoff`. If a dictionary is given, it will apply to the corresponding
            column. If a column cannot be found in the dictionary, it defaults to
            `no`.
            no: Do not filter variants (no binary files are generated in outdir).
            gt: Filter variants with MAF greater than `cutoff`.
            lt: Filter variants with MAF less than `cutoff`.
            ge: Filter variants with MAF greater than or equal to `cutoff`.
            le: Filter variants with MAF less than or equal to `cutoff`.
        plot (flag): If set, plot the distribution of allele frequencies.
        devpars (ns): The device parameters for the plot.
            - width (type=int): Width of the plot
            - height (type=int): Height of the plot
            - res (type=int): Resolution of the plot
    """
    input = "indir:dir"
    output = "outdir:dir:{{in.indir | stem}}.freq"
    lang = config.lang.rscript
    envs = {
        "plink": config.exe.plink2,
        "ncores": config.misc.ncores,
        "modifier": "none",
        "gz": False,
        "cutoff": {},
        "filter": {},
        "plot": True,
        "devpars": {"width": 1000, "height": 800, "res": 100},
    }
    script = "file://../scripts/snp/PlinkFreq.R"
    plugin_opts = {"report": "file://../reports/snp/PlinkFreq.svelte"}


class PlinkUpdateName(Proc):
    """Update variant names in PLINK files.

    See also <https://www.cog-genomics.org/plink/2.0/data#update_map>.

    Input:
        indir: Input directory containing the PLINK files.
            Including `.bed`, `.bim`, and `.fam` files
        namefile: File containing the variant names to update.
            Either a file containing two columns, the first column is the old
            variant name, and the second column is the new variant name.
            Or a VCF file containing the variant names to update.
            When a VCF file is given, the chromosome, position, and reference and
            alternate alleles will be used to match the variants.

    Output:
        outdir: Output directory containing the updated PLINK files.
            Including `.bed`, `.bim`, and `.fam` files

    Envs:
        ncores: Number of cores/threads to use, will pass to plink `--threads` option
        plink: Path to PLINK v2
        bcftools: Path to bcftools
        match_alt (choice): How to match alternate alleles when `in.namefile`
            is a VCF file.
            - exact: Matches alternate alleles exactly.
            - all: Matches alternate alleles regardless of the order.
                `chr1:100:A:T,G` matches `chr1:100:A:G,T` or `chr1:100:A:T,G`.
            - any: Matches any alternate allele.
                For example, `chr1:100:A:T,G` matches `chr1:100:A:G,C`
            - first_included: Matches when the first allele is included.
                For example, `chr1:100:A:T,G` matches `chr1:100:A:C,T`.
            - first: Match first alternate allele
                For example, `chr1:100:A:T,G` matches `chr1:100:A:T`.
            - none: Do not match alternate alleles
    """
    input = "indir:dir, namefile:file"
    output = "outdir:dir:{{in.indir | stem}}.newnames"
    lang = config.lang.python
    envs = {
        "ncores": config.misc.ncores,
        "plink": config.exe.plink2,
        "bcftools": config.exe.bcftools,
        "match_alt": "exact",
    }
    script = "file://../scripts/snp/PlinkUpdateName.py"
