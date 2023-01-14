"""The CNVkit pipeline."""
import pandas
from diot import Diot
from datar.tibble import tibble
from pipen_cli_run import Pipeline

from ..core.config import config


DEFAULT_COLS = Diot(
    group="Group",
    purity="Purity",
    snpvcf="SnpVcf",
    bam="Bam",
    vcf_sample_id="VcfSampleId",
    vcf_normal_id="VcfNormalId",
    sex="Sex",
)


# channel modifier helpers
def _get_metadf(ch):
    return pandas.read_csv(ch.outfile.tolist()[0], sep="\t", header=0)


def _get_all_bams(ch, bamcol):
    return _get_metadf(ch)[bamcol].tolist()


class CNVkitPipeline(Pipeline):
    """The CNVkit pipeline

    Unlike `cnvkit.py batch`, this decouples the steps of the `batch` command so
    that we can control the details of each step.

    The pipeline requires following options:

    Input files:
    - metafile: a tab-separated file (see the next section)
    - baitfile: Potentially targeted genomic regions.
        E.g. all possible exons for the reference genome.
        Format - BED, interval list, etc.
        Optional if `method` is `wgs`.
    - case: The group name of samples in `metacols.group` to call CNVs for.
    - control: The group name of samples in `metacols.group` to use as reference
    - metacols: The column names for each type of information in metafile
        - group: The column name in the metafile that indicates the sample group
        - purity: The column name in the metafile that indicates the sample
            purity
        - snpvcf: The column name in the metafile that indicates the path to
            the SNP VCF file
        - bam: The column name in the metafile that indicates the path to the
            BAM file
        - sample: The column name in the metafile that indicates the sample ID
        - vcf_sample_id: The column name in the metafile that indicates the
            sample ID in the VCF file
        - vcf_normal_id: The column name in the metafile that indicates the
            normal sample ID in the VCF file
        - sex: The column name in the metafile that indicates the sample
            sex

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
        Requires `DNAcopy` to be installed in R, defaults to
        `config.lang.rscript`
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
    - cnvkit: the path to the cnvkit.py executable, defaults to
        `config.exe.cnvkit` from `./.biopipen.toml` or `~/.biopipen.toml`.
    - convert: Linux `convert` command to convert pdf to png
    - access.file: the path to the access file if not given, will use
        options under `access` to generate one.
    - access.exclude: the regions to exclude for generating access file
    - access.min_gap_size: the minimum gap size for generating access file
    - autobin.bp_per_bin: Desired average number of sequencing read bases mapped
        to each bin.
    - autobin.target_max_size: Maximum size of target bins.
    - autobin.target_min_size: Minimum size of target bins.
    - autobin.antitarget_max_size: Maximum size of antitarget bins.
    - autobin.antitarget_min_size: Minimum size of antitarget bins.
    - coverage.min_mapq: Minimum mapping quality for a read to be included.
    - coverage.count: Get read depths by counting read midpoints within each bin
        (An alternative algorithm).
    - reference.cluster: Calculate and store summary stats for clustered subsets
        of the normal samples with similar coverage profiles.
    - reference.min_cluster_size: Minimum cluster size to keep in reference
        profiles
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
    - segment.smooth_cbs: Perform an additional smoothing before CBS
        segmentation, which in some cases may increase the sensitivity.
        Used only for CBS method.
    - scatter.chromosome: Chromosome or chromosomal range, e.g. 'chr1' or
        'chr1:2333000-2444000', to display. If a range is given, all targeted
        genes in this range will be shown, unless -g/--gene is also given.
    - scatter.gene: Name of gene or genes (comma-separated) to display.
    - scatter.range_list: A list of chromosomal ranges to display, as BED,
        interval list or 'chr:start-end' text. Creates focal plots similar to
        -c/--chromosome for each listed region, combined into a multi-page PDF.
    - scatter.width: Width of margin to show around the selected gene(s)
        (-g/--gene) or small chromosomal region (-c/--chromosome).
        [Default: 1000000]
    - scatter.antitarget_marker: Plot antitargets using this symbol when
        plotting in a selected chromosomal region (-g/--gene or -c/--chromosome)
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
    - heatmap_cns.by_bin: Plot data x-coordinates by bin indices instead of
        genomic coordinates. All bins will be shown with equal width,
        no blank regions will be shown, and x-axis values indicate
        bin number (within chromosome) instead of genomic position.
    - heatmap_cns.chromosome: Chromosome (e.g. 'chr1') or chromosomal range
        (e.g. 'chr1:2333000-2444000') to display.
    - heatmap_cns.desaturate: Tweak color saturation to focus on significant
        changes.
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
    - Sample: The sample_id used for target/antitarget files. If not provided,
        the sample_id will be the first part of basename of the bam file.
        For exapmle: `D123.tumor.bam -> D123`
    - `<bam>`: The path to the bam file, better using absolute path.
    - `<group>`: The type of the sample, defining the tumor/normal samples.
    - `<sex>`: Guess each sample from coverage of X and Y chromosomes if
        not given.
    - `<purity>`: Estimated tumor cell fraction, a.k.a. purity or cellularity.
    - `<snpvcf>`: file name containing variants for segmentation by allele
        frequencies.
    - `<vcf_sample_id>`: Sample ID in the VCF file.
    - `<vcf_normal_id>`: Normal sample ID in the VCF file.

    To run this pipeline from command line, with the `pipen-run` plugin:
    >>> # In this case, `pipeline.cnvkit_pipeline.metafile` must be provided
    >>> pipen run cnvkit_pipeline CNVkitPipeline \
    >>>     +config <config.toml> <other pipeline args>

    To use this as a dependency for other pipelines:
    >>> from biopipen.ns.cnvkit_pipeline import build_processes
    >>> MetaFile, CNVkitAccess = build_processes(<options>)
    >>> MetaFile.requires = <the process to generate the metafile>
    """

    defaults = config.pipeline.cnvkit_pipeline

    def col(self, name: str):
        metacols = self.options.get("metacols", {})
        return metacols.get(name, DEFAULT_COLS[name])

    def build(self):
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

        cnvkit = self.options.get("cnvkit", config.exe.cnvkit)
        convert = self.options.get("cnvkit", config.exe.convert)
        ref = self.options.get("reffa", config.ref.reffa)
        rscript = self.options.get("rscript", config.lang.rscript)
        ncores = self.options.get("ncores", config.misc.ncores)

        access_opts = self.options.get("access", {})
        autobin_opts = self.options.get("autobin", {})
        coverage_opts = self.options.get("coverage", {})
        reference_opts = self.options.get("reference", {})
        fix_opts = self.options.get("fix", {})
        segment_opts = self.options.get("segment", {})
        scatter_opts = self.options.get("scatter", {})
        diagram_opts = self.options.get("diagram", {})
        heatmap_cns_opts = self.options.get("heatmap_cns", {})
        heatmap_cnr_opts = self.options.get("heatmap_cnr", {})
        call_opts = self.options.get("call", {})

        class MetaFile(File2Proc):
            """Pass by the metafile"""
            if "metafile" in self.options:
                input_data = [self.options.metafile]

        self.starts.append(MetaFile)
        self.procs.MetaFile = MetaFile

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

        self.starts.append(CNVkitAccess)
        self.procs.CNVkitAccess = CNVkitAccess

        class CNVkitAutobin(CNVkitAutobin):
            requires = [MetaFile, CNVkitAccess]
            input_data = lambda ch1, ch2: [
                (
                    _get_all_bams(ch1, self.col("bam")),
                    ch2.iloc[0, 0],
                    self.options.get("baitfile"),
                ),
            ]
            envs = {
                "cnvkit": cnvkit,
                "method": self.options.get("method", "hybrid"),
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
                "annotate": self.options.get(
                    "annotate", config.ref.refflat
                ),
                "short_names": self.options.get("short_names", False),
                "ref": ref,
            }

        self.procs.CNVkitAutobin = CNVkitAutobin

        class CNVkitCoverageTarget(CNVkitCoverage):
            requires = [MetaFile, CNVkitAutobin]
            input_data = lambda ch1, ch2: tibble(
                _get_all_bams(ch1, self.col("bam")),
                target_file=ch2.target_file.tolist()[0],
            )
            envs = {
                "cnvkit": cnvkit,
                "count": coverage_opts.get("count", False),
                "min_mapq": coverage_opts.get("min_mapq", 0),
                "ncores": ncores,
            }

        self.procs.CNVkitCoverageTarget = CNVkitCoverageTarget

        class CNVkitCoverageAntitarget(CNVkitCoverage):
            requires = [MetaFile, CNVkitAutobin]
            input_data = lambda ch1, ch2: tibble(
                _get_all_bams(ch1, self.col("bam")),
                target_file=ch2.antitarget_file.tolist()[0],
            )
            envs = {
                "cnvkit": cnvkit,
                "count": coverage_opts.get("count", False),
                "min_mapq": coverage_opts.get("min_mapq", 0),
                "ncores": ncores,
            }

        self.procs.CNVkitCoverageAntitarget = CNVkitCoverageAntitarget

        class CNVkitReference(CNVkitReference):
            def _get_input_data(ch1, ch2, ch3, ch4):
                metadf = _get_metadf(ch1)
                normal_masks = metadf[self.col("group")] == self.options.control

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
                        ",".join(metadf[self.col("sex")][normal_masks])
                        if self.col("sex") in metadf.columns
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
                "min_cluster_size": reference_opts.get(
                    "min_cluster_size", False
                ),
                "male_reference": self.options.get("male_reference", False),
                "no_gc": self.options.get("no_gc", False),
                "no_edge": self.options.get("no_edge", False),
                "no_rmask": self.options.get("no_rmask", False),
                "ref": ref,
            }

        self.procs.CNVkitReference = CNVkitReference

        class CNVkitFix(CNVkitFix):
            def _get_input_data(ch1, ch2, ch3, ch4):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case
                return tibble(
                    target_file=ch2.outfile[tumor_masks],
                    antitarget_file=ch3.outfile[tumor_masks],
                    reference=ch4.outfile,
                    sample_id=metadf["Sample"][tumor_masks],
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
                "no_gc": self.options.get("no_gc", False),
                "no_edge": self.options.get("no_edge", False),
                "no_rmask": self.options.get("no_rmask", False),
            }

        self.procs.CNVkitFix = CNVkitFix

        class CNVkitSegment(CNVkitSegment):
            def _get_input_data(ch1, ch2):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case
                return tibble(
                    chrfile=ch2.outfile,
                    vcf=(
                        metadf[self.col("snpvcf")][tumor_masks]
                        if self.col("snpvcf") in metadf.columns
                        else [None]
                    ),
                    sample_id=(
                        metadf[self.col("vcf_sample_id")][tumor_masks]
                        if self.col("vcf_sample_id") in metadf.columns
                        else [None]
                    ),
                    normal_id=(
                        metadf[self.col("vcf_normal_id")][tumor_masks]
                        if self.col("vcf_normal_id") in metadf.columns
                        else [None]
                    ),
                )

            requires = [MetaFile, CNVkitFix]
            input_data = _get_input_data
            envs = {
                "cnvkit": cnvkit,
                "method": segment_opts.get("method", "cbs"),
                "threshold": segment_opts.get("threshold", False),
                "drop_low_coverage": self.options.get(
                    "drop_low_coverage", False
                ),
                "drop_outliers": segment_opts.get("drop_outliers", 10),
                "rscript": rscript,
                "ncores": ncores,
                "smooth_cbs": segment_opts.get("smooth_cbs", False),
                "min_variant_depth": self.options.get("min_variant_depth", 20),
                "zygosity_freq": self.options.get("zygosity_freq", 0.25),
            }

        self.procs.CNVkitSegment = CNVkitSegment

        class CNVkitScatter(CNVkitScatter):
            def _get_input_data(ch1, ch2, ch3):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case

                return tibble(
                    chrfile=ch2.outfile,
                    cnsfile=ch3.outfile,
                    vcf=(
                        metadf[self.col("snpvcf")][tumor_masks]
                        if self.col("snpvcf") in metadf.columns
                        else [None]
                    ),
                    sample_id=(
                        metadf[self.col("vcf_sample_id")][tumor_masks]
                        if self.col("vcf_sample_id") in metadf.columns
                        else [None]
                    ),
                    normal_id=(
                        metadf[self.col("vcf_normal_id")][tumor_masks]
                        if self.col("vcf_normal_id") in metadf.columns
                        else [None]
                    ),
                )

            requires = [MetaFile, CNVkitFix, CNVkitSegment]
            input_data = _get_input_data
            envs = {"cnvkit": cnvkit, "convert": convert, **scatter_opts}

        self.ends.append(CNVkitScatter)
        self.procs.CNVkitScatter = CNVkitScatter

        class CNVkitDiagram(CNVkitDiagram):
            def _get_input_data(ch1, ch2, ch3):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case
                return tibble(
                    chrfile=ch2.outfile,
                    cnsfile=ch3.outfile,
                    sample_sex=(
                        metadf[self.col("sex")][tumor_masks]
                        if self.col("sex") in metadf.columns
                        else [None]
                    ),
                )

            requires = [MetaFile, CNVkitFix, CNVkitSegment]
            input_data = _get_input_data
            envs = {
                "cnvkit": cnvkit,
                "convert": convert,
                **diagram_opts,
            }

        self.ends.append(CNVkitDiagram)
        self.procs.CNVkitDiagram = CNVkitDiagram

        class CNVkitHeatmapCns(CNVkitHeatmap):
            def _get_input_data(ch1, ch2):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case
                return tibble(
                    segfiles=[ch2.outfile.tolist()],
                    sample_sex=(
                        ",".join(metadf[self.col("sex")][tumor_masks])
                        if self.col("sex") in metadf.columns
                        else [None]
                    ),
                )

            requires = [MetaFile, CNVkitSegment]
            input_data = _get_input_data
            envs = {
                "cnvkit": cnvkit,
                "convert": convert,
                "male_reference": self.options.get("male_reference", False),
                **heatmap_cns_opts,
            }

        self.ends.append(CNVkitHeatmapCns)
        self.procs.CNVkitHeatmapCns = CNVkitHeatmapCns

        if heatmap_cnr_opts:

            class CNVkitHeatmapCnr(CNVkitHeatmap):
                def _get_input_data(ch1, ch2):
                    metadf = _get_metadf(ch1)
                    tumor_masks = metadf[self.col("group")] == self.options.case
                    return tibble(
                        segfiles=[ch2.outfile.tolist()],
                        sample_sex=(
                            ",".join(metadf[self.col("sex")][tumor_masks])
                            if self.col("sex") in metadf.columns
                            else [None]
                        ),
                    )

                requires = [MetaFile, CNVkitFix]
                input_data = _get_input_data
                envs = {
                    "cnvkit": cnvkit,
                    "convert": convert,
                    "male_reference": self.options.get(
                        "male_reference", False
                    ),
                    **heatmap_cnr_opts,
                }

            self.ends.append(CNVkitHeatmapCnr)
            self.procs.CNVkitHeatmapCnr = CNVkitHeatmapCnr

        class CNVkitCall(CNVkitCall):
            def _get_input_data(ch1, ch2, ch3):
                metadf = _get_metadf(ch1)
                tumor_masks = metadf[self.col("group")] == self.options.case
                return tibble(
                    cnrfile=ch2.outfile,
                    cnsfile=ch3.outfile,
                    vcf=(
                        metadf[self.col("snpvcf")][tumor_masks]
                        if self.col("snpvcf") in metadf.columns
                        else [None]
                    ),
                    sample_id=(
                        metadf[self.col("vcf_sample_id")][tumor_masks]
                        if self.col("vcf_sample_id") in metadf.columns
                        else [None]
                    ),
                    normal_id=(
                        metadf[self.col("vcf_normal_id")][tumor_masks]
                        if self.col("vcf_normal_id") in metadf.columns
                        else [None]
                    ),
                    sample_sex=(
                        metadf[self.col("sex")][tumor_masks]
                        if self.col("sex") in metadf.columns
                        else [None]
                    ),
                    purity=(
                        metadf[self.col("purity")][tumor_masks]
                        if self.col("purity") in metadf.columns
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
                "thresholds": call_opts.get(
                    "thresholds", "-1.1,-0.25,0.2,0.7"
                ),
                "ploidy": call_opts.get("ploidy", 2),
                "drop_low_coverage": self.options.get(
                    "drop_low_coverage", False
                ),
                "male_reference": self.options.get("male_reference", False),
                "min_variant_depth": self.options.get("min_variant_depth", 20),
                "zygosity_freq": self.options.get("zygosity_freq", 0.25),
            }

        self.ends.append(CNVkitCall)
        self.procs.CNVkitCall = CNVkitCall
