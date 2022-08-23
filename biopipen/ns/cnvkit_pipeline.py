"""The CNVkit pipeline.

Unlike `cnvkit.py batch`, this decouples the steps of the `batch` command so
that we can control the details of each step.

The pipeline requires following options:

Input files:
- metafile: a tab-separated file (see the next section)
- baitfile: Potentially targeted genomic regions.
    E.g. all possible exons for the reference genome.
    Format - BED, interval list, etc.
    Optional if `method` is `wgs`.
- type_col: The column name in the metafile that indicates the sample type.
- type_tumor: The type of tumor samples in `type_col` column of `metafile`
- type_normal: The type of normal samples in `type_col` column of `metafile`

Other command options from `cnvkit.py batch`:
- method: Sequencing protocol: hybridization capture ('hybrid'),
    targeted amplicon sequencing ('amplicon'),
    or whole genome sequencing ('wgs'). Determines
    whether and how to use antitarget bins.
- male_reference: Use or assume a male reference (i.e. female samples
    will have +1 log-CNR of chrX; otherwise male samples would have
    -1 chrX).
- drop_low_coverage: Drop very-low-coverage bins before segmentation to
    avoid false-positive deletions in poor-quality tumor samples.
- ncores: number of cores to use, defaults to `config.misc.ncores`
- rscript: Path to the Rscript excecutable to use for running R code.
    Requires `DNAcopy` to be installed in R, defaults to `config.lang.rscript`
- reffa: the reference genome (e.g. hg19.fa)
- annotate: Use gene models from this file to assign names to the target
    regions. Format: UCSC refFlat.txt or ensFlat.txt file (preferred),
    or BED, interval list, GFF, or similar.
- short_names: Reduce multi-accession bait labels to be short and
    consistent.
- refcnn: Reuse the reference if created before. Otherwise use options under
    `reference` to create a new reference.
- no_gc: Skip GC correction for `cnvkit.py reference/fix`.
- no_edge: Skip edge-effect correction for `cnvkit.py reference/fix`.
- no_rmask: Skip RepeatMasker correction for `cnvkit.py reference/fix`.
- min_variant_depth: Minimum read depth for a SNV to be displayed
    in the b-allele frequency plot.
- zygosity_freq: Ignore VCF's genotypes (GT field) and instead infer
    zygosity from allele frequencies.
- cnvkit: the path to the cnvkit.py executable, defaults to `config.exe.cnvkit`
    from `./.biopipen.toml` or `~/.biopipen.toml`.
- convert: Linux `convert` command to convert pdf to png
- access.file: the path to the access file if not given, will use
    options under `access` to generate one.
- access.exclude: the regions to exclude for generating access file
- access.min_gap_size: the minimum gap size for generating access file
- autobin.bp_per_bin: Desired average number of sequencing read bases mapped to
    each bin.
- autobin.target_max_size: Maximum size of target bins.
- autobin.target_min_size: Minimum size of target bins.
- autobin.antitarget_max_size: Maximum size of antitarget bins.
- autobin.antitarget_min_size: Minimum size of antitarget bins.
- coverage.min_mapq: Minimum mapping quality for a read to be included.
- coverage.count: Get read depths by counting read midpoints within each bin.
    (An alternative algorithm).
- reference.cluster: Calculate and store summary stats for clustered subsets of
    the normal samples with similar coverage profiles.
- reference.min_cluster_size: Minimum cluster size to keep in reference profiles
- fix.sample_id: True to use `Sample` from metafile
- fix.cluster: Compare and use cluster-specific values present in the
    reference profile. Requires reference.cluster=True.
- segment.method: Method to use for segmentation.
    Candidates - cbs, flasso, haar, none, hmm, hmm-tumor, hmm-germline
- segment.threshold: Significance threshold
    (p-value or FDR, depending on method)
    to accept breakpoints during segmentation. For HMM methods,
    this is the smoothing window size.
- segment.drop_outliers: Drop outlier bins more than this many multiples of
    the 95th quantile away from the average within a rolling window.
    Set to 0 for no outlier filtering.
- segment.smooth_cbs: Perform an additional smoothing before CBS segmentation,
    which in some cases may increase the sensitivity. Used only for CBS method.
- scatter.chromosome: Chromosome or chromosomal range, e.g. 'chr1' or
    'chr1:2333000-2444000', to display. If a range is given, all targeted genes
    in this range will be shown, unless -g/--gene is also given.
- scatter.gene: Name of gene or genes (comma-separated) to display.
- scatter.range_list: A list of chromosomal ranges to display, as BED,
    interval list or 'chr:start-end' text. Creates focal plots similar to
    -c/--chromosome for each listed region, combined into a multi-page PDF.
- scatter.width: Width of margin to show around the selected gene(s) (-g/--gene)
    or small chromosomal region (-c/--chromosome). [Default: 1000000]
- scatter.antitarget_marker: Plot antitargets using this symbol when plotting
    in a selected chromosomal region (-g/--gene or -c/--chromosome).
    Same as targets if not given
- scatter.by_bin: Plot data x-coordinates by bin indices instead of genomic
    coordinates. All bins will be shown with equal width, no blank
    regions will be shown, and x-axis values indicate bin number
    (within chromosome) instead of genomic position.
- scatter.segment_color: Plot segment lines in this color. Value can be
    any string accepted by matplotlib, e.g. 'red' or '#CC0000'.
- scatter.trend: Draw a smoothed local trendline on the scatter plot.
- scatter.y_max: y-axis upper limit.
- scatter.y_min: y-axis lower limit.
- scatter.title: Plot title. Sample ID if not provided.
- scatter.convert_args: The default arguments for convert
- scatter.cases: If more than one case to plot
- diagram.threshold: Copy number change threshold to label genes.
- diagram.min_probes: Minimum number of covered probes to label a gene.
- diagram.no_shift_xy: Don't adjust the X and Y chromosomes according
    to sample sex.
- diagram.title: Plot title. Sample ID if not provided.
- diagram.convert_args: The default arguments for convert
- diagram.cases: If more than one case to plot
- heatmap_cns.by_bin: Plot data x-coordinates by bin indices instead of genomic
    coordinates. All bins will be shown with equal width,
    no blank regions will be shown, and x-axis values indicate
    bin number (within chromosome) instead of genomic position.
- heatmap_cns.chromosome: Chromosome (e.g. 'chr1') or chromosomal range
    (e.g. 'chr1:2333000-2444000') to display.
- heatmap_cns.desaturate: Tweak color saturation to focus on significant changes.
- heatmap_cns.no_shift_xy: Don't adjust the X and Y chromosomes according to
    sample sex.
- heatmap_cns.convert_args: The default arguments for convert
- heatmap_cns.cases: If more than one case to plot
- heatmap_chr: Similar as `heatmap_cns`, but use the `.cnr` file
    (from `cnvkit.py fix`)
- call.center: Re-center the log2 ratio values using this estimator of
    the center or average value.
- call.center_at: Subtract a constant number from all log2 ratios.
    For "manual" re-centering, in case the --center option gives
    unsatisfactory results.)
- call.filter: Merge segments flagged by the specified filter(s) with
    the adjacent segment(s).
- call.method: Calling method (threshold, clonal or none).
- call.thresholds: Hard thresholds for calling each integer copy number,
    separated by commas.
- call.ploidy: Ploidy of the sample cells.
- call.purity: Estimated tumor cell fraction, a.k.a. purity or cellularity.
-
A metafile should be with the following columns:
- Sample: The sample_id used for target/antitarget files. If not provided, the
    sample_id will be the first part of basename of the bam file.
    For exapmle: `D123.tumor.bam -> D123`
- BamFile: The path to the bam file, better using absolute path.
- `<type_col>`: The type of the sample, defining the tumor/normal samples.
- SampleSex: Guess each sample from coverage of X and Y chromosomes if
    not given.
- VcfFile: file name containing variants for segmentation by allele frequencies.
- VcfSampleId: Sample ID in the VCF file.
- VcfNormalId: Normal sample ID in the VCF file.

To run this pipeline from command line, with the `pipen-run` plugin:
>>> # In this case, `pipeline.cnvkit_pipeline.metafile` must be provided
>>> pipen run cnvkit_pipeline main +config <config.toml> <other pipeline args>

To use this as a dependency for other pipelines:
>>> from biopipen.ns.cnvkit_pipeline import build_processes
>>> MetaFile, CNVkitAccess = build_processes(<options>)
>>> MetaFile.requires = <the process to generate the metafile>
"""
from typing import Mapping, Any

import pandas
from diot import Diot
from pipen import Proc, Pipen
from datar.tibble import tibble

DEFAULT = {

}


def build_processes(options: Mapping[str, Any] = None) -> Proc:
    """Build processes for CNVkit pipeline"""
    from ..core.config import config
    from .cnvkit import (
        CNVkitAccess,
        CNVkitAutobin,
        CNVkitCoverage,
        CNVkitReference,
        CNVkitFix,
        CNVkitSegment,
        CNVkitScatter,
        CNVkitDiagram,
        CNVkitHeatmap,
        CNVkitCall,
    )
    from .misc import File2Proc

    options = (
        Diot(DEFAULT)
        | (
            config
            .get("pipeline", {})
            .get("cnvkit_pipeline", {})
        )
        | (options or {})
    )

    cnvkit = options.get("cnvkit", config.exe.cnvkit)
    convert = options.get("cnvkit", config.exe.convert)
    ref = options.get("reffa", config.ref.reffa)
    rscript = options.get("rscript", config.lang.rscript)
    ncores = options.get("ncores", config.misc.ncores)

    # channel modifier helpers
    def _get_metadf(ch):
        return pandas.read_csv(ch.outfile.tolist()[0], sep="\t", header=0)

    def _get_all_bams(ch):
        return _get_metadf(ch)["BamFile"].tolist()

    class MetaFile(File2Proc):
        """Pass by the metafile"""
        if "metafile" in options:
            input_data = [options.metafile]

    access_opts = options.get("access", {})
    if "file" not in access_opts:
        class CNVkitAccess(CNVkitAccess):
            input_data = [access_opts.get("exclude", [])]
            envs = {
                "cnvkit": cnvkit,
                "min_gap_size": access_opts.get("min_gap_size", 5000),
                "ref": ref,
            }
    else:
        class CNVkitAccess(File2Proc):
            input_data = [access_opts.file]

    autobin_opts = options.get("autobin", {})
    class CNVkitAutobin(CNVkitAutobin):
        requires = [MetaFile, CNVkitAccess]
        input_data = lambda ch1, ch2: [
            (_get_all_bams(ch1), ch2.iloc[0, 0], options.get("baitfile")),
        ]
        envs = {
            "cnvkit": cnvkit,
            "method": options.get("method", "hybrid"),
            "bp_per_bin": autobin_opts.get("bp_per_bin", 100000),
            "target_max_size": autobin_opts.get("target_max_size", 20000),
            "target_min_size": autobin_opts.get("target_min_size", 20),
            "antitarget_max_size": autobin_opts.get(
                "antitarget_max_size",
                500000,
            ),
            "antitarget_min_size": autobin_opts.get(
                "antitarget_min_size",
                1000,
            ),
            "annotate": options.get("annotate", config.ref.refflat),
            "short_names": options.get("short_names", False),
            "ref": ref,
        }

    coverage_opts = options.get("coverage", {})
    class CNVkitCoverageTarget(CNVkitCoverage):
        requires = [MetaFile, CNVkitAutobin]
        input_data = lambda ch1, ch2: tibble(
            _get_all_bams(ch1),
            target_file=ch2.target_file.tolist()[0],
        )
        envs = {
            "cnvkit": cnvkit,
            "count": coverage_opts.get("count", False),
            "min_mapq": coverage_opts.get("min_mapq", 0),
            "ncores": ncores,
        }

    class CNVkitCoverageAntitarget(CNVkitCoverage):
        requires = [MetaFile, CNVkitAutobin]
        input_data = lambda ch1, ch2: tibble(
            _get_all_bams(ch1),
            target_file=ch2.antitarget_file.tolist()[0],
        )
        envs = {
            "cnvkit": cnvkit,
            "count": coverage_opts.get("count", False),
            "min_mapq": coverage_opts.get("min_mapq", 0),
            "ncores": ncores,
        }

    reference_opts = options.get("reference", {})
    class CNVkitReference(CNVkitReference):
        def _get_input_data(ch1, ch2, ch3, ch4):
            metadf = _get_metadf(ch1)
            normal_masks = metadf[options.type_col] == options.type_normal
            return tibble(
                covfiles=(
                    [None]
                    if sum(normal_masks) == 0
                    else [
                        ch2.outfile[normal_masks].tolist()
                        + ch3.outfile[normal_masks].tolist()
                    ]
                ),
                target_file=ch4.target_file,
                antitarget_file=ch4.antitarget_file,
                sample_sex=(
                    ",".join(metadf.SampleSex[normal_masks])
                    if "SampleSex" in metadf.columns
                    else [None]
                ),
            )

        requires = [
            MetaFile,
            CNVkitCoverageTarget,
            CNVkitCoverageAntitarget,
            CNVkitAutobin,
        ]
        input_data = _get_input_data
        envs = {
            "cnvkit": cnvkit,
            "cluster": reference_opts.get("cluster", False),
            "min_cluster_size": reference_opts.get("min_cluster_size", False),
            "male_reference": options.get("male_reference", False),
            "no_gc": options.get("no_gc", False),
            "no_edge": options.get("no_edge", False),
            "no_rmask": options.get("no_rmask", False),
            "ref": ref,
        }

    fix_opts = options.get("fix", {})
    class CNVkitFix(CNVkitFix):
        def _get_input_data(ch1, ch2, ch3, ch4):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor
            return tibble(
                target_file=ch2.outfile[tumor_masks],
                antitarget_file=ch3.outfile[tumor_masks],
                reference=ch4.outfile,
                sample_id=(
                    metadf.SampleID[tumor_masks]
                    if "SampleID" in metadf.columns
                    else [None]
                ),
            )

        requires = [
            MetaFile,
            CNVkitCoverageTarget,
            CNVkitCoverageAntitarget,
            CNVkitReference,
        ]
        input_data = _get_input_data
        envs = {
            "cnvkit": cnvkit,
            "cluster": fix_opts.get("cluster", False),
            "no_gc": options.get("no_gc", False),
            "no_edge": options.get("no_edge", False),
            "no_rmask": options.get("no_rmask", False),
        }

    segment_opts = options.get("segment", {})
    class CNVkitSegment(CNVkitSegment):
        def _get_input_data(ch1, ch2):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor
            return tibble(
                chrfile=ch2.outfile,
                vcf=(
                    metadf.SnpVcf[tumor_masks]
                    if "SnpVcf" in metadf.columns
                    else [None]
                ),
                sample_id=(
                    metadf.VcfSampleID[tumor_masks]
                    if "VcfSampleID" in metadf.columns
                    else [None]
                ),
                normal_id=(
                    metadf.NormalID[tumor_masks]
                    if "NormalID" in metadf.columns
                    else [None]
                ),
            )

        requires = [MetaFile, CNVkitFix]
        input_data = _get_input_data
        envs = {
            "cnvkit": cnvkit,
            "method": segment_opts.get("method", "cbs"),
            "threshold": segment_opts.get("threshold", False),
            "drop_low_coverage": options.get("drop_low_coverage", False),
            "drop_outliers": segment_opts.get("drop_outliers", 10),
            "rscript": rscript,
            "ncores": ncores,
            "smooth_cbs": segment_opts.get("smooth_cbs", False),
            "min_variant_depth": options.get("min_variant_depth", 20),
            "zygosity_freq": options.get("zygosity_freq", 0.25),
        }

    scatter_opts = options.get("scatter", {})
    class CNVkitScatter(CNVkitScatter):
        def _get_input_data(ch1, ch2, ch3):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor

            return tibble(
                chrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                vcf=(
                    metadf.SnpVcf[tumor_masks]
                    if "SnpVcf" in metadf.columns
                    else [None]
                ),
                sample_id=(
                    metadf.VcfSampleID[tumor_masks]
                    if "VcfSampleID" in metadf.columns
                    else [None]
                ),
                normal_id=(
                    metadf.NormalID[tumor_masks]
                    if "NormalID" in metadf.columns
                    else [None]
                ),
            )

        requires = [MetaFile, CNVkitFix, CNVkitSegment]
        input_data = _get_input_data
        envs={
            "cnvkit": cnvkit,
            "convert": convert,
            **scatter_opts
        }

    diagram_opts = options.get("diagram", {})
    class CNVkitDiagram(CNVkitDiagram):
        def _get_input_data(ch1, ch2, ch3):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor
            return tibble(
                chrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                sample_sex=(
                    metadf.SampleSex[tumor_masks]
                    if "SampleSex" in metadf.columns
                    else [None]
                ),
            )

        requires = [MetaFile, CNVkitFix, CNVkitSegment]
        input_data = _get_input_data
        envs={
            "cnvkit": cnvkit,
            "convert": convert,
            **diagram_opts,
        }

    heatmap_cns_opts = options.get("heatmap_cns", {})
    class CNVkitHeatmapCns(CNVkitHeatmap):
        def _get_input_data(ch1, ch2):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor
            return tibble(
                segfiles=[ch2.outfile[tumor_masks].tolist()],
                sample_sex=(
                    ",".join(metadf.SampleSex[tumor_masks])
                    if "SampleSex" in metadf.columns
                    else [None]
                ),
            )

        requires = [MetaFile,  CNVkitSegment]
        input_data = _get_input_data
        envs={
            "cnvkit": cnvkit,
            "convert": convert,
            "male_reference": options.get("male_reference", False),
            **heatmap_cns_opts,
        }

    heatmap_cnr_opts = options.get("heatmap_cnr", {})
    if heatmap_cnr_opts:
        class CNVkitHeatmapCnr(CNVkitHeatmap):
            def _get_input_data(ch1, ch2):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[options.type_col] == options.type_tumor
                return tibble(
                    segfiles=[ch2.outfile[tumor_masks].tolist()],
                    sample_sex=(
                        ",".join(metadf.SampleSex[tumor_masks])
                        if "SampleSex" in metadf.columns
                        else [None]
                    ),
                )

            requires = [MetaFile,  CNVkitFix]
            input_data = _get_input_data
            envs={
                "cnvkit": cnvkit,
                "convert": convert,
                "male_reference": options.get("male_reference", False),
                **heatmap_cnr_opts,
            }

    call_opts = options.get("call", {})
    class CNVkitCall(CNVkitCall):
        def _get_input_data(ch1, ch2, ch3):
            metadf = _get_metadf(ch1)
            tumor_masks = metadf[options.type_col] == options.type_tumor
            return tibble(
                cnrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                vcf=(
                    metadf.SnpVcf[tumor_masks]
                    if "SnpVcf" in metadf.columns
                    else [None]
                ),
                sample_id=(
                    metadf.VcfSampleID[tumor_masks]
                    if "VcfSampleID" in metadf.columns
                    else [None]
                ),
                normal_id=(
                    metadf.VcfNormalID[tumor_masks]
                    if "NormalID" in metadf.columns
                    else [None]
                ),
                sample_sex=(
                    metadf.SampleSex[tumor_masks]
                    if "SampleSex" in _get_metadf(ch1).columns
                    else [None]
                ),
            )
        requires = [MetaFile, CNVkitFix, CNVkitSegment]
        input_data = _get_input_data
        envs = {
            "cnvkit": cnvkit,
            "center": call_opts.get("center", "median"),
            "center_at": call_opts.get("center_at", False),
            "filter": call_opts.get("filter", False),
            "method": call_opts.get("method", "threshold"),
            "thresholds": call_opts.get("thresholds", "-1.1,-0.25,0.2,0.7"),
            "ploidy": call_opts.get("ploidy", 2),
            "purity": call_opts.get("purity", False),
            "drop_low_coverage": options.get("drop_low_coverage", False),
            "male_reference": options.get("male_reference", False),
            "min_variant_depth": options.get("min_variant_depth", 20),
            "zygosity_freq": options.get("zygosity_freq", 0.25),
        }

    return MetaFile, CNVkitAccess

def main() -> Pipen:
    """Build a pipeline for `pipen run` to run"""
    return Pipen(
        name="cnvkit-pipeline",
        desc="CNV calling pipeline using cnvkit",
    ).set_start(build_processes())
