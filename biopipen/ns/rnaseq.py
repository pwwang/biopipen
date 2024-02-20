"""RNA-seq data analysis"""

from ..core.proc import Proc
from ..core.config import config


class UnitConversion(Proc):
    """Convert expression value units back and forth

    See <https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>
    and <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm>.

    Following converstions are supported -
    * `count -> cpm, fpkm/rpkm, fpkmuq/rpkmrq, tpm, tmm`
    * `fpkm/rpkm -> count, tpm, cpm`
    * `tpm -> count, fpkm/rpkm, cpm`
    * `cpm -> count, fpkm/rpkm, tpm`
    NOTE that during some conversions, `sum(counts/effLen)` is approximated to
    `sum(counts)/sum(effLen) * length(effLen))`

    You can also use this process to just transform the expression values, e.g., take
    log2 of the expression values. In this case, you can set `inunit` and `outunit` to
    `count` and `log2(count + 1)` respectively.

    Input:
        infile: Input file containing expression values
            The file should be a matrix with rows representing genes and columns
            representing samples.
            It could be an RDS file containing a data frame or a matrix, or a
            text file containing a matrix with tab as the delimiter. The text
            file can be gzipped.

    Output:
        outfile: Output file containing the converted expression values
            The file will be a matrix with rows representing genes and columns
            representing samples.

    Envs:
        inunit: The input unit of the expression values.
            You can also use an expression to indicate the input unit, e.g.,
            `log2(counts + 1)`. The expression should be like `A * fn(B*X + C) + D`,
            where `A`, `B`, `C` and `D` are constants, `fn` is a function, and X is
            the input unit.
            Currently only `expr`, `sqrt`, `log2`, `log10` and `log` are supported as
            functions.
            Supported input units are:
            * counts/count/rawcounts/rawcount: raw counts.
            * cpm: counts per million.
            * fpkm/rpkm: fragments per kilobase of transcript per million.
            * fpkmuq/rpkmuq: upper quartile normalized FPKM/RPKM.
            * tpm: transcripts per million.
            * tmm: trimmed mean of M-values.
        outunit: The output unit of the expression values. An expression can also be
            used for transformation (e.g. `log2(tpm + 1)`). If `inunit` is `count`,
            then this means we are converting raw counts to tpm, and transforming it
            to `log2(tpm + 1)` as the output. Any expression supported by `R` can be
            used. Same units as `inunit` are supported.
        refexon: Path to the reference exon gff file.
        meanfl (type=auto): A file containing the mean fragment length for each sample
            by rows (samples as rowname), without header.
            Or a fixed universal estimated number (1 used by TCGA).
        nreads (type=auto): The estimatied total number of reads for each sample.
            or you can pass a file with the number for each sample by rows
            (samples as rowname), without header.
            When converting `fpkm/rpkm -> count`, it should be total reads of that sample.
            When converting `cpm -> count`: it should be total reads of that sample.
            When converting `tpm -> count`: it should be total reads of that sample.
            When converting `tpm -> cpm`: it should be total reads of that sample.
            When converting `tpm -> fpkm/rpkm`: it should be `sum(fpkm)` of that sample.
            It is not used when converting `count -> cpm, fpkm/rpkm, tpm`.
    """  # noqa: E501
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.rscript
    envs = {
        "inunit": None,
        "outunit": None,
        "refexon": config.ref.refexon,
        "meanfl": 1,
        "nreads": 1_000_000,
    }
    script = "file://../scripts/rnaseq/UnitConversion.R"


class Simulation(Proc):
    """Simulate RNA-seq data using ESCO/RUVcorr package

    Input:
        ngenes: Number of genes to simulate
        nsamples: Number of samples to simulate
            If you want to force the process to re-simulate for the same
            `ngenes` and `nsamples`, you can set a different value for `envs.seed`.
            Note that the samples will be shown as cells in the output (since
            the simulation is designed for single-cell RNA-seq data).

    Output:
        outfile: Output file containing the simulated data with rows representing
            genes and columns representing samples.
        outdir: Output directory containing the simulated data
            `sim.rds` and `True.rds` will be generated.
            For `ESCO`, `sim.rds` contains the simulated data in a
            `SingleCellExperiment` object, and `True.rds` contains the matrix of true
            counts.
            For `RUVcorr`, `sim.rds` contains the simulated data in list with
            `Truth`, A matrix containing the values of Xβ; `Y` A matrix containing the
            values in `Y`; `Noise` A matrix containing the values in `Wα`; `Sigma`
            A matrix containing the true gene-gene correlations, as defined by Xβ; and
            `Info` A matrix containing some of the general information about the
            simulation.
            For all matrices, rows represent genes and columns represent samples.

    Envs:
        tool (choice): Which tool to use for simulation.
            - ESCO: uses the [ESCO](https://github.com/JINJINT/ESCO) package.
            - RUVcorr: uses the [RUVcorr](https://rdrr.io/bioc/RUVcorr/) package.
        ncores (type=int): Number of cores to use.
        seed (type=int): Random seed.
            If not set, seed will not be set.
        esco_args (ns): Additional arguments to pass to the simulation function.
            - save (choice): Which type of data to save to `out.outfile`.
                - `simulated-truth`: saves the simulated true counts.
                - `zero-inflated`: saves the zero-inflated counts.
                - `down-sampled`: saves the down-sampled counts.
            - type (choice): Which type of heterogenounity to use.
                - single: produces a single population.
                - group: produces distinct groups.
                - tree: produces distinct groups but admits a tree structure.
                - traj: produces distinct groups but admits a smooth trajectory
                    structure.
            - <more>: See <https://rdrr.io/github/JINJINT/ESCO/man/escoParams.html>.
        ruvcorr_args (ns): Additional arguments to pass to the simulation
            function.
            - <more>: See <https://rdrr.io/bioc/RUVcorr/man/simulateGEdata.html>.
        transpose_output (flag): If set, the output will be transposed.
        index_start (type=int): The index to start from when naming the samples.
            Affects the sample names in `out.outfile` only.
    """
    input = "ngenes:var, nsamples:var"
    output = [
        "outfile:file:{{in.ngenes}}x{{in.nsamples}}.sim/simulated.txt",
        "outdir:dir:{{in.ngenes}}x{{in.nsamples}}.sim",
    ]
    lang = config.lang.rscript
    envs = {
        "tool": "RUVcorr",
        "ncores": config.misc.ncores,
        "type": "single",
        "esco_args": {
            "dropout-type": "none",
            "save": "simulated-truth",
            "type": "single",
        },
        "ruvcorr_args": {},
        "seed": None,
        "transpose_output": False,
        "index_start": 1,
    }
    script = "file://../scripts/rnaseq/Simulation.R"
