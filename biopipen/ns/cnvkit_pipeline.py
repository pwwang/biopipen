"""The CNVkit pipeline."""
from __future__ import annotations
from typing import TYPE_CHECKING, Any
from functools import lru_cache

import pandas
from diot import Diot
from datar.tibble import tibble
from biopipen.core.proc import Proc
from pipen_cli_run import Pipeline, process

from ..core.config import config

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property

if TYPE_CHECKING:
    from pandas import DataFrame

DEFAULT_COLS = Diot(
    group="Group",
    purity="Purity",
    snpvcf="SnpVcf",
    bam="Bam",
    vcf_sample_id="VcfSampleId",
    vcf_normal_id="VcfNormalId",
    sex="Sex",
    guess_baits="GuessBaits",
)

DEFAULT_OPTS = Diot(
    metafile=None,
    baitfile=None,
    accfile=None,
    cnvkit=config.exe.cnvkit,
    convert=config.exe.convert,
    rscript=config.lang.rscript,
    samtools=config.exe.samtools,
    ncores=config.misc.ncores,
    reffa=config.ref.reffa,
    annotate=config.ref.refflat,
    short_names=True,
    method="hybrid",
    guessbaits=False,
    heatmap_cnr=False,
    case=None,
    control=None,
    access_excludes=[],
    guessbaits_guided=None,
    male_reference=False,
    drop_low_coverage=False,
    min_variant_depth=20,
    no_gc=False,
    no_edge=False,
    no_rmask=False,
    zygosity_freq=0.25,
)


@lru_cache()
def _metadf(metafile: str) -> DataFrame:
    return pandas.read_csv(metafile, sep="\t", header=0)


def _1st(df: DataFrame) -> Any:
    return df.iloc[0, 0]


class _MetaCol:
    """Get the column name from the metafile"""
    def __init__(self, cols, default_cols):
        self.cols = cols or {}
        self.default_cols = default_cols

    def __getattr__(self, name):
        return self.cols.get(name, self.default_cols[name])


class CNVkitPipeline(Pipeline):
    """The CNVkit pipeline

    Unlike `cnvkit.py batch`, this decouples the steps of the `batch` command so
    that we can control the details of each step.

    The pipeline requires following options:

    Input files:
    - metafile: a tab-separated file (see the next section)
    - baitfile: Potentially targeted genomic regions.
        E.g. all possible exons for the reference genome.
        This is optional when `method` is `wgs`.
    - accfile: The accessible genomic regions.
        If not given, use `cnvkit.py access` to generate one.

    Special options:
    - access_excludes: File(s) with regions to be excluded for
        `cnvkit.py access`.
    - guessbaits_guided: Whether to use guided mode for guessing baits.
    - guessbaits: Guess the bait file from the bam files, either guided or
        unguided.
        If False, `baitfile` is used. Otherwise, if `baitfile` is given, use it
        (guided), otherwise use `accfile` (unguided).
        The bam files with `metacols.guess_baits` column set to `True`, `TRUE`,
        `true`, `1`, `Yes`, `YES`, or `yes` will be used to guess the bait file.
    - heatmap_cnr: Whether to generate a heatmap of the .cnr files
        (bin-level signals). This is allowed to set to False, it will take
        longer to run.
    - case: The group name of samples in `metacols.group` to call CNVs for.
        If not specified, use all samples. In such a case, `control` must not be
        specified, as we are using a flat reference.
    - control: The group name of samples in `metacols.group` to use as reference
        if not specified, use a flat reference.
    - metacols: The column names for each type of information in metafile
        - group: The column name in the metafile that indicates the sample group
            Default: `Group`
        - purity: The column name in the metafile that indicates the sample
            purity. Default: `Purity`
        - snpvcf: The column name in the metafile that indicates the path to
            the SNP VCF file. Default: `SnpVcf`
        - bam: The column name in the metafile that indicates the path to the
            BAM file. Default: `Bam`
        - vcf_sample_id: The column name in the metafile that indicates the
            sample ID in the VCF file. Default: `VcfSampleId`
        - vcf_normal_id: The column name in the metafile that indicates the
            normal sample ID in the VCF file. Default: `VcfNormalId`
        - sex: The column name in the metafile that indicates the sample
            sex. Default: `Sex`
        - guess_baits: The column name in the metafile that indicates whether
            to guess the bait file from the bam files. Default: `GuessBaits`

    Global options that are used by multiple processes
    (can be overriden individually by `[<proc>.envs.xxx]`):
    - cnvkit: the path to the cnvkit.py executable, defaults to
        `config.exe.cnvkit` from `./.biopipen.toml` or `~/.biopipen.toml`.
    - rscript: Path to the Rscript excecutable to use for running R code.
        Requires `DNAcopy` to be installed in R, defaults to
        `config.lang.rscript`
    - samtools: Path to samtools, used for guessing bait file.
    - convert: Linux `convert` command to convert pdf to png
        So that they can be embedded in the HTML report.
    - ncores: number of cores to use, defaults to `config.misc.ncores`
    - reffa: the reference genome (e.g. hg19.fa)
        Used by `CNVkitAccess`, `CNVkitAutobin` and `CNVkitReference`
    - annotate: Use gene models from this file to assign names to the
        target regions. Format: UCSC refFlat.txt or ensFlat.txt file
        (preferred), or BED, interval list, GFF, or similar.
    - short_names: Reduce multi-accession bait labels to be short and consistent
    - method: Sequencing protocol: hybridization capture ('hybrid'),
        targeted amplicon sequencing ('amplicon'),
        or whole genome sequencing ('wgs'). Determines
        whether and how to use antitarget bins.
    - male_reference: Use or assume a male reference (i.e. female samples
        will have +1 log-CNR of chrX; otherwise male samples would have
        -1 chrX).
        Used by `CNVkitReference`, `CNVkitCall`, `CNVkitHeatmapCns` and
        `CNVkitHeatmapCnr`.
    - drop_low_coverage: Drop very-low-coverage bins before segmentation to
        avoid false-positive deletions in poor-quality tumor samples.
        Used by `CNVkitSegment` and `CNVkitCall`
    - no_gc: Skip GC correction for `cnvkit.py reference/fix`.
    - no_edge: Skip edge-effect correction for `cnvkit.py reference/fix`.
    - no_rmask: Skip RepeatMasker correction for `cnvkit.py reference/fix`.
        no_* options are used by `CNVkitReference` and `CNVkitFix`
    - min_variant_depth: Minimum read depth for a SNV to be displayed
        in the b-allele frequency plot.
        Used by `CNVkitSegment` and `CNVkitCall`
    - zygosity_freq: Ignore VCF's genotypes (GT field) and instead infer
        zygosity from allele frequencies.
        Used by `CNVkitSegment` and `CNVkitCall`

    Options for different processes can be specified by `[CNVkitXXX.envs.xxx]`
    See `biopipen.ns.cnvkit.CNVkitXXX` for more details.

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
    - `<guess_baits>`: Guess the bait file from the bam file

    To run this pipeline from command line, with the `pipen-run` plugin:
    >>> # In this case, `pipeline.cnvkit_pipeline.metafile` must be provided
    >>> pipen run cnvkit_pipeline CNVkitPipeline \
    >>>     +config <config.toml> <other pipeline args>

    To use this as a dependency for other pipelines:
    >>> from biopipen.ns.cnvkit_pipeline import CNVkitPipeline
    >>> pipeline = CNVkitPipeline(<options>)
    >>> # pipeline.starts: Start processes of the pipeline
    >>> # pipeline.ends: End processes of the pipeline
    >>> # pipeline.procs.<proc>: The process with name <proc>
    """

    defaults = config.pipeline.cnvkit_pipeline

    @cached_property
    def col(self):
        """Get the column names by self.col.<colname>"""
        return _MetaCol(self.options.get("metacols"), DEFAULT_COLS)

    @process(start=True)
    def build_metafile(self):
        """Build MetaFile process"""
        from .misc import File2Proc

        class MetaFile(File2Proc):
            """Pass by the metafile"""
            # Do not require metafile, as we could use the pipeline as part of
            # another pipeline, which can generate a metafile
            # Remember to set the dependency in the pipeline:
            # >>> pipeline.procs.MetaFile.requires = [other_pipeline.procs]
            # where other_pipeline.procs generate the metafile
            if self.options.metafile:
                input_data = [self.options.metafile]

        return MetaFile

    @process(start=True)
    def build_cnvkit_access(self):
        """Build CNVkitAccess process"""
        if self.options.get("accfile"):
            from .misc import File2Proc

            class CNVkitAccess(File2Proc):
                input_data = [self.options.accfile]
        else:
            from .cnvkit import CNVkitAccess

            excludes = self.options.get("excludes", [])
            if not isinstance(excludes, (list, tuple)):
                excludes = [excludes]

            class CNVkitAccess(CNVkitAccess):
                # can be overwritten by [CNVkitAccess.in.exludes]
                input_data = [excludes]
                envs = {
                    "cnvkit": self.options.cnvkit,
                    "ref": self.options.reffa,
                }

        return CNVkitAccess

    @process
    def build_cnvkit_guessbaits(self, metafile_proc, access_proc):
        """Build CNVkitGuessBaits process"""
        from .cnvkit import CNVkitGuessBaits

        if self.options.guessbaits_guided is None:
            raise ValueError(
                "`guessbaits.guided` must be specified, expecting True or False"
            )

        def _guess_baits_bams(ch):
            df = _metadf(_1st(ch))
            if self.col.guess_baits not in df:
                # Use all bams
                return df.loc[:, self.col.bam].tolist()

            # Use only specified
            guess_baits = df[self.col.guess_baits]
            return df.loc[
                guess_baits == True  # noqa
                | guess_baits == "True"
                | guess_baits == "TRUE"
                | guess_baits == "true"
                | guess_baits == "1"
                | guess_baits == 1
                | guess_baits == "yes"
                | guess_baits == "YES"
                | guess_baits == "Yes",
                self.col.bam,
            ].tolist()

        if self.options.guessbaits_guided:
            if not self.options.baitfile:
                raise ValueError(
                    "`baitfile` must be specified for guided mode "
                    "to guess baits. See: "
                    "https://cnvkit.readthedocs.io/en/stable/scripts.html"
                )

            class CNVkitGuessBaits(CNVkitGuessBaits):
                requires = metafile_proc
                input_data = lambda metafile_ch: tibble(
                    bamfiles=[_guess_baits_bams(metafile_ch)],
                    atfile=self.options.baitfile,
                )
                envs = {
                    "cnvkit": self.options.cnvkit,
                    "samtools": self.options.samtools,
                    "ncores": self.options.ncores,
                    "ref": self.options.reffa,
                    "guided": True,
                }
        else:  # unguided
            class CNVkitGuessBaits(CNVkitGuessBaits):
                requires = metafile_proc, access_proc
                input_data = lambda metafile_ch, access_ch: tibble(
                    bamfiles=[_guess_baits_bams(metafile_ch)],
                    accessfile=_1st(access_ch),
                )
                envs = {
                    "cnvkit": self.options.cnvkit,
                    "samtools": self.options.samtools,
                    "ncores": self.options.ncores,
                    "ref": self.options.reffa,
                    "guided": False,
                }

        return CNVkitGuessBaits

    @process
    def build_cnvkit_autobin(self, metafile_proc, access_proc, guessbaits_proc):
        """Build CNVkitAutobin process"""
        from .cnvkit import CNVkitAutobin

        class CNVkitAutobin(CNVkitAutobin):
            if guessbaits_proc:
                requires = metafile_proc, access_proc, guessbaits_proc
                input_data = lambda ch1, ch2, ch3: tibble(
                    bamfiles=[_metadf(_1st(ch1))[self.col.bam].tolist()],
                    accfile=_1st(ch2),
                    baitfile=(
                        _1st(ch3)
                        if self.options.guessbaits
                        else self.options.baitfile
                    ),
                )
            else:
                requires = metafile_proc, access_proc
                input_data = lambda ch1, ch2: tibble(
                    bamfiles=[_metadf(_1st(ch1))[self.col.bam].tolist()],
                    accfile=_1st(ch2),
                    baitfile=self.options.baitfile,
                )
            envs = {
                "cnvkit": self.options.cnvkit,
                "method": self.options.method,
                "annotate": self.options.annotate,
                "short_names": self.options.short_names,
                "ref": self.options.reffa,
            }

        return CNVkitAutobin

    @process
    def build_cnvkit_coverage(self, metafile_proc, autobin_proc, anti):
        """Build CNVkitTargetCoverage and CNVkitAntiTargetCoverage processes"""
        from .cnvkit import CNVkitCoverage

        return Proc.from_proc(
            CNVkitCoverage,
            name="CNVkitCoverageAnittarget" if anti else "CNVkitCoverageTarget",
            requires=[metafile_proc, autobin_proc],
            input_data=lambda ch1, ch2: tibble(
                _metadf(_1st(ch1))[self.col.bam].tolist(),
                target_file=ch2[
                    "antitarget_file" if anti else "target_file"
                ].tolist()[0],
            ),
            envs={
                "cnvkit": self.options.cnvkit,
                "ncores": self.options.ncores,
                "ref": self.options.reffa,
            }
        )

    @process
    def build_cnvkit_reference(
        self,
        metafile_proc,
        target_coverage_proc,
        antitarget_coverage_proc,
        autobin_proc,
    ):
        """Build CNVkitReference process"""
        from .cnvkit import CNVkitReference

        def _input_data(ch1, ch2, ch3, ch4):
            metadf = _metadf(_1st(ch1))

            if self.options.control:
                # Use control samples to build reference
                control_masks = metadf[self.col.group] == self.options.control
                covfiles = [
                    ch2.outfile[control_masks].tolist()
                    + ch3.outfile[control_masks].tolist()
                ]
                target_file = None
                antitarget_file = None
                if self.col.sex in metadf:
                    sample_sex = ",".join(metadf[self.col.sex][control_masks])
                else:
                    sample_sex = [None]
            else:
                # Build a flat reference
                covfiles = [None]
                target_file = ch4.target_file
                antitarget_file = ch4.antitarget_file
                sample_sex = [None]

            return tibble(
                covfiles=covfiles,
                target_file=target_file,
                antitarget_file=antitarget_file,
                sample_sex=sample_sex,
            )

        class CNVkitReference(CNVkitReference):
            requires = [
                metafile_proc,
                target_coverage_proc,
                antitarget_coverage_proc,
                autobin_proc,
            ]
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "no_gc": self.options.no_gc,
                "no_edge": self.options.no_edge,
                "no_rmask": self.options.no_rmask,
                "ref": self.options.reffa,
            }

        return CNVkitReference

    @process
    def build_cnvkit_fix(
        self,
        metafile_proc,
        target_coverage_proc,
        antitarget_coverage_proc,
        reference_proc,
    ):
        """Build CNVkitFix process"""
        from .cnvkit import CNVkitFix

        if not self.options.case and self.options.control:
            raise ValueError(
                "`case` is not specified, meaning using all samples as cases, "
                "but `control` is specified (we can only use a flat reference "
                "in this case)."
            )

        def _input_data(ch1, ch2, ch3, ch4):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                target_file=ch2.outfile[tumor_masks],
                antitarget_file=ch3.outfile[tumor_masks],
                reference=ch4.outfile,
                sample_id=metadf["Sample"][tumor_masks],
            )

        class CNVkitFix(CNVkitFix):
            requires = [
                metafile_proc,
                target_coverage_proc,
                antitarget_coverage_proc,
                reference_proc,
            ]
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "no_gc": self.options.no_gc,
                "no_edge": self.options.no_edge,
                "no_rmask": self.options.no_rmask,
            }

        return CNVkitFix

    @process
    def build_cnvkit_segment(self, metafile_proc, fix_proc):
        """Build CNVkitSegment process"""
        from .cnvkit import CNVkitSegment

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                chrfile=ch2.outfile,
                vcf=(
                    metadf[self.col.snpvcf][tumor_masks]
                    if self.col.snpvcf in metadf
                    else [None]
                ),
                sample_id=(
                    metadf[self.col.vcf_sample_id][tumor_masks]
                    if self.col.vcf_sample_id in metadf
                    else [None]
                ),
                normal_id=(
                    metadf[self.col.vcf_normal_id][tumor_masks]
                    if self.col.vcf_normal_id in metadf.columns
                    else [None]
                ),
            )

        class CNVkitSegment(CNVkitSegment):
            requires = metafile_proc, fix_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "rscript": self.options.rscript,
                "ncores": self.options.ncores,
            }

        return CNVkitSegment

    @process(end=True)
    def build_cnvkit_scatter(self, metafile_proc, fix_proc, segment_proc):
        """Build CNVkitScatter process"""
        from .cnvkit import CNVkitScatter

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                chrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                vcf=(
                    metadf[self.col.snpvcf][tumor_masks]
                    if self.col.snpvcf in metadf
                    else [None]
                ),
                sample_id=(
                    metadf[self.col.vcf_sample_id][tumor_masks]
                    if self.col.vcf_sample_id in metadf
                    else [None]
                ),
                normal_id=(
                    metadf[self.col.vcf_normal_id][tumor_masks]
                    if self.col.vcf_normal_id in metadf
                    else [None]
                ),
            )

        class CNVkitScatter(CNVkitScatter):
            requires = metafile_proc, fix_proc, segment_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "convert": self.options.convert,
            }

        return CNVkitScatter

    @process(end=True)
    def build_cnvkit_diagram(self, metafile_proc, fix_proc, segment_proc):
        """Build CNVkitDiagram process"""
        from .cnvkit import CNVkitDiagram

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                chrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                sample_sex=(
                    metadf[self.col.sex][tumor_masks]
                    if self.col.sex in metadf
                    else [None]
                ),
            )

        class CNVkitDiagram(CNVkitDiagram):
            requires = metafile_proc, fix_proc, segment_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "convert": self.options.convert,
            }

        return CNVkitDiagram

    @process(end=True)
    def build_cnvkit_heatmap_cns(self, metafile_proc, segment_proc):
        """Build CNVkitHeatmapCns process"""
        from .cnvkit import CNVkitHeatmap

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                segfiles=[ch2.outfile.tolist()],
                sample_sex=(
                    ",".join(metadf[self.col.sex][tumor_masks])
                    if self.col.sex in metadf
                    else [None]
                ),
            )

        class CNVkitHeatmapCns(CNVkitHeatmap):
            """Heatmap of segment-level signals of multiple samples"""
            requires = metafile_proc, segment_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "convert": self.options.convert,
                "male_reference": self.options.male_reference,
            }

        return CNVkitHeatmapCns

    @process(end=True)
    def build_cnvkit_heatmap_cnr(self, metafile_proc, fix_proc):
        """Build CNVkitHeatmapCnr process"""
        from .cnvkit import CNVkitHeatmap

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                segfiles=[ch2.outfile.tolist()],
                sample_sex=(
                    ",".join(metadf[self.col.sex][tumor_masks])
                    if self.col.sex in metadf
                    else [None]
                ),
            )

        class CNVkitHeatmapCnr(CNVkitHeatmap):
            """Heatmap of bin-level signals of multiple samples"""
            requires = metafile_proc, fix_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "convert": self.options.convert,
                "male_reference": self.options.male_reference,
            }

        return CNVkitHeatmapCnr

    @process(end=True)
    def build_cnvkit_call(self, metafile_proc, fix_proc, segment_proc):
        """Build CNVkitCall process"""
        from .cnvkit import CNVkitCall

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.options.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.options.case

            return tibble(
                cnrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                vcf=(
                    metadf[self.col.snpvcf][tumor_masks]
                    if self.col.snpvcf in metadf
                    else [None]
                ),
                sample_id=(
                    metadf[self.col.vcf_sample_id][tumor_masks]
                    if self.col.vcf_sample_id in metadf
                    else [None]
                ),
                normal_id=(
                    metadf[self.col.vcf_normal_id][tumor_masks]
                    if self.col.vcf_normal_id in metadf
                    else [None]
                ),
                sample_sex=(
                    metadf[self.col.sex][tumor_masks]
                    if self.col.sex in metadf
                    else [None]
                ),
                purity=(
                    metadf[self.col.purity][tumor_masks]
                    if self.col.purity in metadf
                    else [None]
                ),
            )

        class CNVkitCall(CNVkitCall):
            requires = metafile_proc, fix_proc, segment_proc
            input_data = _input_data
            envs = {
                "cnvkit": self.options.cnvkit,
                "drop_low_coverage": self.options.drop_low_coverage,
                "male_reference": self.options.male_reference,
                "min_variant_depth": self.options.min_variant_depth,
                "zygosity_freq": self.options.zygosity_freq,
            }

        return CNVkitCall

    def build(self):
        self.options = DEFAULT_OPTS | self.options

        MetaFile = self.build_metafile()
        CNVkitAccess = self.build_cnvkit_access()

        CNVkitGuessBaits = None
        if self.options.guessbaits:
            CNVkitGuessBaits = self.build_cnvkit_guessbaits(
                MetaFile,
                CNVkitAccess,
            )

        CNVkitAutobin = self.build_cnvkit_autobin(
            MetaFile,
            CNVkitAccess,
            CNVkitGuessBaits,
        )

        CNVkitCoverageTarget = self.build_cnvkit_coverage(
            MetaFile,
            CNVkitAutobin,
            anti=False,
        )
        CNVkitCoverageAntitarget = self.build_cnvkit_coverage(
            MetaFile,
            CNVkitAutobin,
            anti=True,
        )
        CNVkitReference = self.build_cnvkit_reference(
            MetaFile,
            CNVkitCoverageTarget,
            CNVkitCoverageAntitarget,
            CNVkitAutobin,
        )
        CNVkitFix = self.build_cnvkit_fix(
            MetaFile,
            CNVkitCoverageTarget,
            CNVkitCoverageAntitarget,
            CNVkitReference,
        )
        CNVkitSegment = self.build_cnvkit_segment(MetaFile, CNVkitFix)

        # end processes
        self.build_cnvkit_scatter(MetaFile, CNVkitFix, CNVkitSegment)
        self.build_cnvkit_diagram(MetaFile, CNVkitFix, CNVkitSegment)
        self.build_cnvkit_heatmap_cns(MetaFile, CNVkitSegment)
        if self.options.heatmap_cnr:
            self.build_cnvkit_heatmap_cnr(MetaFile, CNVkitFix)
        self.build_cnvkit_call(MetaFile, CNVkitFix, CNVkitSegment)
