"""CNV/CNA-related processes, mostly tertiary analysis"""

from ..core.proc import Proc
from ..core.config import config


class AneuploidyScore(Proc):
    """Chromosomal arm SCNA/aneuploidy

    The CAAs in this process are calculated using Cohen-Sharir method
    See https://github.com/quevedor2/aneuploidy_score

    Input:
        segfile: The seg file, generally including chrom, start, end and
            seg.mean (the log2 ratio)

    Output:
        outdir: The output directory containing the CAAs, AS and a histogram
            plot to show the CAAs for each chromosome arm

    Requires:
        - name: AneuploidyScore
          check: |
            {{proc.lang}} <(echo "library(AneuploidyScore)")
        - name: ucsc.hg19.cytoband
          if: {{ proc.envs.genome == 'hg19' }}
          check: |
            {{proc.lang}} <(echo "library(ucsc.hg19.cytoband)")
        - name: ucsc.hg38.cytoband
          if: {{ proc.envs.genome == 'hg38' }}
          check: |
            {{proc.lang}} <(echo "library(ucsc.hg38.cytoband)")
    """
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
        "include_sex": False,
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
        group_col: The column name in the metafile to group the samples
        heatmap_cases: The cases to be included in the heatmap
            By default, all arms are included. If specified, keys are the names
            of the cases and values are the arms, which will be included in
            the heatmap. The list of arms should be a subset of `chr<N>_p` and
            `chr<N>_q`, where `<N>` is the chromosome number from 1 to 22, X, Y.
            You can also use `ALL` to include all arms.
    """
    input = "asdirs:dirs, metafile:file"
    output = (
        "outdir:dir:{{in.asdirs | first | stem}}_etc.aneuploidy_score_summary"
    )
    lang = config.lang.rscript
    script = "file://../scripts/cnv/AneuploidyScoreSummary.R"
    envs = {"group_col": None, "heatmap_cases": {"All-Arms": "ALL"}}
    plugin_opts = {
        "report": "file://../reports/cnv/AneuploidyScoreSummary.svelte",
    }
