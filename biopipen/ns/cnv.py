"""CNV/CNA-related processes, mostly tertiary analysis"""

from ..core.proc import Proc
from ..core.config import config


class AneuploidyScore(Proc):
    """Chromosomal arm SCNA/aneuploidy

    The CAAs in this process are calculated using Cohen-Sharir method
    See https://github.com/quevedor2/aneuploidy_score

    Input:
        segfile: The seg file, generally including chrom, start, end and
            seg.mean (the log2 ratio).
            It is typically a tab-delimited file or a BED file.
            If so, envs.chrom_col, envs.start_col, envs.end_col and envs.seg_col
            are the 1st, 2nd, 3rd and 5th columns, respectively.
            It can also be a VCF file. If so, envs.chrom_col and envs.start_col
            are not required.
            `end_col` and `envs.seg_col` will be a field in the INFO column.
            [`VariantAnnotation`](https://rdrr.io/bioc/VariantAnnotation/)
            is required to extract the INFO field.

    Output:
        outdir: The output directory containing the CAAs, AS and a histogram
            plot to show the CAAs for each chromosome arm

    Envs:
        chrom_col: The column name for chromosome
        start_col: The column name for start position
        end_col: The column name for end position
        seg_col: The column name for seg.mean
        cn_col: The column name for copy number
        segmean_transform (text): A R function to transform `seg.mean`
            The transformed value will be used to calculate the CAAs
        cn_transform (type=auto): A R function to transform `seg.mean` into
            copy number, or a list of cutoffs to determine the copy number.
            See https://cnvkit.readthedocs.io/en/stable/pipeline.html#calling-methods.
            If this is give, `cn_col` will be ignored.
        genome: The genome version, hg19 or hg38
        threshold (type=float): The threshold to determine whether a chromosome
            arm is gained or lost.
        wgd_gf (type=float): The fraction of the genome that is affected by WGD
        excl_chroms (list): The chromosomes to be excluded
            Works with/without `chr` prefix.

    Requires:
        AneuploidyScore:
            - check: {{proc.lang}} <(echo "library(AneuploidyScore)")
        ucsc.hg19.cytoband:
            - if: {{ proc.envs.genome == 'hg19' }}
            - check: {{proc.lang}} <(echo "library(ucsc.hg19.cytoband)")
        ucsc.hg38.cytoband:
            - if: {{ proc.envs.genome == 'hg38' }}
            - check: {{proc.lang}} <(echo "library(ucsc.hg38.cytoband)")
    """  # noqa: E501
    input = "segfile:file"
    output = "outdir:dir:{{in.segfile | stem}}.aneuploidy_score"
    lang = config.lang.rscript
    envs = {
        "chrom_col": "chrom",
        "start_col": "loc.start",
        "end_col": "loc.end",
        "seg_col": "seg.mean",
        "cn_col": None,
        "segmean_transform": None,
        "cn_transform": None,
        "genome": config.ref.genome,
        "threshold": 0.5,
        "wgd_gf": 0.5,
        "excl_chroms": ['chrX', 'chrY'],
    }
    script = "file://../scripts/cnv/AneuploidyScore.R"
    plugin_opts = {
        "report": "file://../reports/cnv/AneuploidyScore.svelte",
        "report_paging": 10,
    }


class AneuploidyScoreSummary(Proc):
    """Summary table and plots from AneuploidyScore

    Input:
        asdirs: The output directories from AneuploidyScore
        metafile: The metafile containing the sample information

    Output:
        outdir: The output directory containing the summary table and plots

    Envs:
        group_cols (type=auto): The column name in the metafile to group the
            samples.
            We also support multiple columns, e.g. `["group1", "group2"]`
            You can also use `group1,group2` to add a secondary grouping
            based on `group2` within each `group1` (only works for 2 groups)
        heatmap_cases (type=json): The cases to be included in the heatmap
            By default, all arms are included. If specified, keys are the names
            of the cases and values are the arms, which will be included in
            the heatmap. The list of arms should be a subset of `chr<N>_p` and
            `chr<N>_q`, where `<N>` is the chromosome number from 1 to 22, X, Y.
            You can also use `ALL` to include all arms.
        sample_name (text): An R function to extract the sample name from
            the file stem (not including `.aneuploidy_score` part)
    """
    input = "asdirs:dirs, metafile:file"
    output = (
        "outdir:dir:{{in.asdirs | first | stem}}_etc.aneuploidy_score_summary"
    )
    lang = config.lang.rscript
    script = "file://../scripts/cnv/AneuploidyScoreSummary.R"
    envs = {
        "group_cols": None,
        "heatmap_cases": {"All-Arms": "ALL"},
        "sample_name": None,
    }
    plugin_opts = {
        "report": "file://../reports/cnv/AneuploidyScoreSummary.svelte",
    }


class TMADScore(Proc):
    """Trimmed Median Absolute Deviation (TMAD) score for CNV

    Reference:
        Mouliere, Chandrananda, Piskorz and Moore et al. Enhanced detection of
        circulating tumor DNA by fragment size analysis Science Translational
        Medicine (2018).

    Input:
        segfile: The seg file, two columns are required:
            * chrom: The chromosome name, used for filtering
            * seg.mean: The log2 ratio.
            It is typically a tab-delimited file or a BED file.
            If so, envs.chrom_col and envs.seg_col
            are the 1st and 5th columns, respectively.
            It can also be a VCF file. If so, envs.chrom_col and envs.start_col
            are not required.
            `end_col` and `envs.seg_col` will be a field in the INFO column.
            [`VariantAnnotation`](https://rdrr.io/bioc/VariantAnnotation/)
            is required to extract the INFO field.

    Output:
        outfile: The output file containing the TMAD score

    Envs:
        chrom_col: The column name for chromosome
        seg_col: The column name for seg.mean
        segmean_transform: The transformation function for seg.mean
        excl_chroms (list): The chromosomes to be excluded
    """
    input = "segfile:file"
    output = "outfile:file:{{in.segfile | stem}}.tmad.txt"
    lang = config.lang.rscript
    envs = {
        "chrom_col": "chrom",
        "seg_col": "seg.mean",
        "segmean_transform": None,
        "excl_chroms": ["chrX", "chrY"],
    }
    script = "file://../scripts/cnv/TMADScore.R"


class TMADScoreSummary(Proc):
    """Summary table and plots for TMADScore

    Input:
        tmadfiles: The output files from TMADScore
        metafile: The metafile containing the sample information
            The first column must be the sample ID

    Output:
        outdir: The output directory containing the summary table and plots

    Envs:
        group_cols (type=auto): The column name in the metafile to group the
            samples Could also be a list of column names
            If not specified, samples will be plotted individually as a barplot
            We also support multiple columns, e.g. `["group1", "group2"]`
            You can also use `group1,group2` to add a secondary grouping
            based on `group2` within each `group1` (only works for 2 groups)
        sample_name (text): An R function to extract the sample name from
            the file stem (not including `.tmad.txt` part)
    """
    input = "tmadfiles:files, metafile:file"
    output = "outdir:dir:{{in.tmadfiles | first | stem0}}_etc.tmad_summary"
    lang = config.lang.rscript
    script = "file://../scripts/cnv/TMADScoreSummary.R"
    envs = {"group_cols": None, "sample_name": None}
    plugin_opts = {
        "report": "file://../reports/cnv/TMADScoreSummary.svelte",
    }
