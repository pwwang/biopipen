"""Provides processes for the regulatory related"""

from ..core.proc import Proc
from ..core.config import config


class MotifScan(Proc):
    """Scan the input sequences for binding sites using motifs.

    Currently only [fimo](https://meme-suite.org/meme/tools/fimo) from MEME suite
    is supported, based on the research/comparisons done by the following reference.

    Reference:
        - [Evaluating tools for transcription factor binding site prediction](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6889335/)

    Input:
        motiffile: File containing motif names.
            The file contains the motif and regulator names.
            The motif names should match the names in the motif database.
            This file must have a header.
            If multiple columns are present, it should be delimited by tab.
        seqfile: File containing sequences in FASTA format.

    Output:
        outdir: Directory containing the results.
            Especially `fimo_output.txt` extending from `fimo.tsv`, which contains:
            1. the results with the regulator information if `envs.regulator_col`
                is provided, otherwise, the `regulator` columns will be filled with
                the motif names.
            2. the original sequence from the fasta file (in.seqfile)
            3. corrected genomic coordinates if the genomic coordinates are included
                in the sequence names.

            See also the `Output` section of
            <https://meme-suite.org/meme/doc/fimo.html>.
            Note that `--no-pgc` is passed to fimo to not parse the genomic coordinates
            from the sequence names by fimo. When fimo parses the genomic coordinates,
            `DDX11L1` in `>DDX11L1::chr1:11869-14412` will be lost.
            The purpose of this is to keep the sequence names as they are in the output.
            If the sequence names are in the format of `>NAME::chr1:START-END`, we will
            correct the coordinates in the output.
            Also note that it requires meme/fimo v5.5.5+ to do this
            (where the --no-pgc option is available).

    Envs:
        tool (choice): The tool to use for scanning.
            Currently only fimo is supported.
            - fimo: Use fimo from MEME suite.
        fimo: The path to fimo binary.
        motif_col: The column name in the motif file containing the motif names.
        regulator_col: The column name in the motif file containing the regulator names.
            Both `motif_col` and `regulator_col` should be the direct column names or
            the index (1-based) of the columns.
            If no `regulator_col` is provided, no regulator information is written in
            the output.
        notfound (choice): What to do if a motif is not found in the database.
            - error: Report error and stop the process.
            - ignore: Ignore the motif and continue.
        motifdb: The path to the motif database. This is required.
            It should be in the format of MEME motif database.
            Databases can be downloaded here: <https://meme-suite.org/meme/doc/download.html>.
            See also introduction to the databases: <https://meme-suite.org/meme/db/motifs>.
        cutoff (type=float): The cutoff for p-value to write the results.
            When `envs.q_cutoff` is set, this is applied to the q-value.
            This is passed to `--thresh` in fimo.
        q (flag): Calculate q-value.
            When `False`, `--no-qvalue` is passed to fimo.
            The q-value calculation is that of Benjamini and Hochberg (BH) (1995).
        q_cutoff (flag): Apply `envs.cutoff` to q-value.
        args (ns): Additional arguments to pass to the tool.
            - <more>: Additional arguments for fimo.
                See: <https://meme-suite.org/meme/doc/fimo.html>
    """  # noqa: E501
    input = "motiffile:file, seqfile:file"
    output = "outdir:dir:{{in.motiffile | stem}}.fimo"
    lang = config.lang.python
    envs = {
        "tool": "fimo",
        "fimo": config.exe.fimo,
        "motif_col": 1,
        "regulator_col": None,
        "notfound": "error",
        "motifdb": config.tf_motifdb,
        "cutoff": 1e-4,
        "q": False,
        "q_cutoff": False,
        "args": {},
    }
    script = "file://../scripts/regulatory/MotifScan.py"


class MotifAffinityTest(Proc):
    """Test the affinity of motifs to the sequences and the affinity change
    due the mutations.

    See also <https://simon-coetzee.github.io/motifBreakR> and
    <https://www.bioconductor.org/packages/release/bioc/vignettes/atSNP/inst/doc/atsnp-vignette.html>

    When using atSNP, motifBreakR is also required to plot the variants and motifs.

    Input:
        motiffile: File containing motif names.
            The file contains the motif and regulator names.
            The motif names should match the names in the motif database.
            This file must have a header.
            If multiple columns are present, it should be delimited by tab.
        varfile: File containing the variants.
            It could be a VCF file or a BED-like file.
            If it is a VCF file, it does not need to be indexed. Only records with `PASS` in the `FILTER` column are used.
            If it is a BED-like file, it should contain the following columns, `chrom`, `start`, `end`, `name`, `score`, `strand`, `ref`, `alt`.

    Output:
        outdir: Directory containing the results.
            For motifBreakR, `motifbreakr.txt` will be created. Records with effect `strong`/`weak` are written (`neutral` is not).
            For atSNP, `atsnp.txt` will be created. Records with p-value (`envs.atsnp_args.p`) < `envs.cutoff` are written.

    Envs:
        ncores (type=int): The number of cores to use.
        tool (choice): The tool to use for the test.
            - motifbreakr: Use motifBreakR.
            - motifBreakR: Use motifBreakR.
            - atsnp: Use atSNP.
            - atSNP: Use atSNP.
        bcftools: The path to bcftools binary.
            Used to convert the VCF file to the BED file when the input is a VCF file.
        motif_col: The column name in the motif file containing the motif names.
            If this is not provided, `envs.regulator_col` and `envs.regmotifs` are required,
            which are used to infer the motif names from the regulator names.
        regulator_col: The column name in the motif file containing the regulator names.
            Both `motif_col` and `regulator_col` should be the direct column names or
            the index (1-based) of the columns.
            If no `regulator_col` is provided, no regulator information is written in
            the output. Otherwise, the regulator information is written in the output in
            the `Regulator` column.
        var_col: The column names in the `in.motiffile` containing the variant information.
            It has to be matching the names in the `in.varfile`. This is helpful when
            we only need to test the pairs of variants and motifs in the `in.motiffile`.
        notfound (choice): What to do if a motif is not found in the database,
            or a regulator is not found in the regulator-motif mapping (envs.regmotifs)
            file.
            - error: Report error and stop the process.
            - ignore: Ignore the motif and continue.
        motifdb: The path to the motif database. This is required.
            It should be in the format of MEME motif database.
            Databases can be downloaded here: <https://meme-suite.org/meme/doc/download.html>.
            See also introduction to the databases: <https://meme-suite.org/meme/db/motifs>.
            [universalmotif](https://github.com/bjmt/universalmotif) is required to read the motif database.
        genome: The genome assembly.
            Used to fetch the sequences around the variants by package, for example, `BSgenome.Hsapiens.UCSC.hg19` is required if
            `hg19`. If it is an organism other than human, please specify the full name of the package, for example, `BSgenome.Mmusculus.UCSC.mm10`.
        cutoff (type=float): The cutoff for p-value to write the results.
        devpars (ns): The default device parameters for the plot.
            - width (type=int): The width of the plot.
            - height (type=int): The height of the plot.
            - res (type=int): The resolution of the plot.
        plot_nvars (type=int): Number of variants to plot.
            Plot top `<plot_nvars>` variants with the largest `abs(alleleDiff)` (motifBreakR) or smallest p-values (atSNP).
        plots (type=json): Specify the details for the plots.
            When specified, `plot_nvars` is ignored.
            The keys are the variant names and the values are the details for the plots, including:
            devpars: The device parameters for the plot to override the default (envs.devpars).
            which: An expression passed to `subset(results, subset = ...)` to get the motifs for the variant to plot.
                Or an integer to get the top `which` motifs.
                For example, `effect == "strong"` to get the motifs with strong effect in motifBreakR result.
        regmotifs: The path to the regulator-motif mapping file.
            It must have header and the columns `Motif` or `Model` for motif names and
            `TF`, `Regulator` or `Transcription factor`  for regulator names.
        motifbreakr_args (ns): Additional arguments to pass to motifBreakR.
            - method (choice): The method to use.
                See details of <https://rdrr.io/bioc/motifbreakR/man/motifbreakR.html>
                and <https://simon-coetzee.github.io/motifBreakR/#methods>.
                - default: Use the default method.
                - log: Use the standard summation of log probabilities
                - ic: Use information content
                - notrans: Use the default method without transformation
        atsnp_args (ns): Additional arguments to pass to atSNP.
            - padj_cutoff (flag): The `envs.cutoff` will be applied to the adjusted p-value.
                Only works for `atSNP`.
            - padj (choice): The method to adjust the p-values.
                Only works for `atSNP`
                - holm: Holm's method
                - hochberg: Hochberg's method
                - hommel: Hommel's method
                - bonferroni: Bonferroni method
                - BH: Benjamini & Hochberg's method
                - BY: Benjamini & Yekutieli's method
                - fdr: False discovery rate
                - none: No adjustment
            - p (choice): Which p-value to use for adjustment and cutoff.
                - pval_ref: p-value for the reference allele affinity score.
                - pval_snp: p-value for the SNP allele affinity score.
                - pval_cond_ref: and
                - pval_cond_snp: conditional p-values for the affinity scores of the reference and SNP alleles.
                - pval_diff: p-value for the affinity score change between the two alleles.
                - pval_rank: p-value for the rank test between the two alleles.
    """  # noqa: E501
    input = "motiffile:file, varfile:file"
    output = "outdir:dir:{{in.motiffile | stem}}.{{envs.tool | lower}}"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "tool": "atsnp",
        "bcftools": config.exe.bcftools,
        "motif_col": None,
        "regulator_col": None,
        "var_col": None,
        "notfound": "error",
        "motifdb": config.ref.tf_motifdb,
        "regmotifs": config.ref.tf_motifs,
        "genome": config.ref.genome,
        "cutoff": 0.05,
        "devpars": {"width": None, "height": None, "res": 100},
        "plot_nvars": 10,
        "plots": {},
        "motifbreakr_args": {"method": "default"},
        "atsnp_args": {"padj_cutoff": True, "padj": "BH", "p": "pval_diff"},
    }
    script = "file://../scripts/regulatory/MotifAffinityTest.R"


class VariantMotifPlot(Proc):
    """A plot with a genomic region surrounding a genomic variant, and
    potentially disrupted motifs.

    Currently only SNVs are supported.

    Input:
        infile: File containing the variants and motifs.
            It is a TAB-delimited file with the following columns:
            - chrom: The chromosome of the SNV. Alias: chr, seqnames.
            - start: The start position of the SNV, no matter 0- or 1-based.
            - end: The end position of the SNV, which will be used as the position of the SNV.
            - strand: Indicating the direction of the surrounding sequence matching the motif.
            - SNP_id: The name of the SNV.
            - REF: The reference allele of the SNV.
            - ALT: The alternative allele of the SNV.
            - providerId: The motif id. It can be specified by `envs.motif_col`.
            - providerName: The name of the motif provider. Optional.
            - Regulator: The regulator name. Optional, can be specified by `envs.regulator_col`.
            - motifPos: The position of the motif, relative to the position of the SNV.
                For example, '-8, 4' means the motif is 8 bp upstream and 4 bp downstream of the SNV.

    Envs:
        genome: The genome assembly.
            Used to fetch the sequences around the variants by package, for example, `BSgenome.Hsapiens.UCSC.hg19` is required if
            `hg19`. If it is an organism other than human, please specify the full name of the package, for example, `BSgenome.Mmusculus.UCSC.mm10`.
        motifdb: The path to the motif database. This is required.
            It should be in the format of MEME motif database.
            Databases can be downloaded here: <https://meme-suite.org/meme/doc/download.html>.
            See also introduction to the databases: <https://meme-suite.org/meme/db/motifs>.
            [universalmotif](https://github.com/bjmt/universalmotif) is required to read the motif database.
        motif_col: The column name in the motif file containing the motif names.
            If this is not provided, `envs.regulator_col` and `envs.regmotifs` are required,
            which are used to infer the motif names from the regulator names.
        regulator_col: The column name in the motif file containing the regulator names.
            Both `motif_col` and `regulator_col` should be the direct column names or
            the index (1-based) of the columns.
            If no `regulator_col` is provided, no regulator information is written in
            the output. Otherwise, the regulator information is written in the output in
            the `Regulator` column.
        regmotifs: The path to the regulator-motif mapping file.
            It must have header and the columns `Motif` or `Model` for motif names and
            `TF`, `Regulator` or `Transcription factor`  for regulator names.
        notfound (choice): What to do if a motif is not found in the database,
            or a regulator is not found in the regulator-motif mapping (envs.regmotifs)
            file.
            - error: Report error and stop the process.
            - ignore: Ignore the motif and continue.
        devpars (ns): The default device parameters for the plot.
            - width (type=int): The width of the plot.
            - height (type=int): The height of the plot.
            - res (type=int): The resolution of the plot.
        plot_vars (type=auto): The variants (SNP_id) to plot.
            A list of variant names to plot or a string with the variant names separated by comma.
            When not specified, all variants are plotted.
    """  # noqa: E501
    input = "infile:file"
    output = "outdir:dir:{{in.infile | stem}}.vmplots"
    lang = config.lang.rscript
    envs = {
        "genome": config.ref.genome,
        "motifdb": config.ref.tf_motifdb,
        "motif_col": "providerId",
        "regulator_col": None,
        "regmotifs": config.ref.tf_motifs,
        "notfound": "error",
        "devpars": {"width": 800, "height": None, "res": 100},
        "plot_vars": None,
    }
    script = "file://../scripts/regulatory/VariantMotifPlot.R"
