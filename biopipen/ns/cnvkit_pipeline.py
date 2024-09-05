"""The CNVkit pipeline."""
from __future__ import annotations
from typing import TYPE_CHECKING, Any
from functools import lru_cache

import pandas
from diot import Diot
from datar.tibble import tibble
from pipen.utils import mark, is_loading_pipeline
from biopipen.core.proc import Proc
from pipen_annotate import annotate
from pipen_args.procgroup import ProcGroup

from ..core.config import config

from functools import cached_property

if TYPE_CHECKING:
    from pandas import DataFrame


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


class CNVkitPipeline(ProcGroup):
    """The CNVkit pipeline

    Unlike `cnvkit.py batch`, this decouples the steps of the `batch` command so
    that we can control the details of each step.

    Options for different processes can be specified by `[CNVkitXXX.envs.xxx]`
    See `biopipen.ns.cnvkit.CNVkitXXX` for more details.

    To run this pipeline from command line, with the `pipen-run` plugin:
    >>> # In this case, `pipeline.cnvkit_pipeline.metafile` must be provided
    >>> pipen run cnvkit_pipeline CNVkitPipeline <other pipeline args>

    To use this as a dependency for other pipelines -
    >>> from biopipen.ns.cnvkit_pipeline import CNVkitPipeline
    >>> pipeline = CNVkitPipeline(<options>)
    >>> # pipeline.starts: Start processes of the pipeline
    >>> # pipeline.ends: End processes of the pipeline
    >>> # pipeline.procs.<proc>: The process with name <proc>

    See also the docs for details
    <https://pwwang.github.io/biopipen/pipelines/cnvkit_pipeline/>

    Args:
        metafile (order=-99): A tab-separated file.
            * Sample: Unique IDs of the samples. Required.
            * `<bam>`: The path to the bam file, better using absolute path.
            * `<group>`: The type of the sample, defining the tumor/normal
                samples.
            * `<sex>`: Guess each sample from coverage of X and Y chromosomes
                if not given.
            * `<purity>`: Estimated tumor cell fraction, a.k.a. purity or
                cellularity.
            * `<snpvcf>`: file name containing variants for segmentation by
                allele frequencies.
            * `<vcf_sample_id>`: Sample ID in the VCF file.
            * `<vcf_normal_id>`: Normal sample ID in the VCF file.
            * `<guess_baits>`: Whether use this bam file to guess the baits
        metacols (ns;order=-98): The column names for each type of information
            in metafile.
            - group (default=Group): The column name in the metafile that
                indicates the sample group
            - purity: The column name in the metafile that indicates the sample
                purity
            - snpvcf: The column name in the metafile that indicates the path to
                the SNP VCFflag
            - bam: The column name in the metafile that indicates the path to
                the BAM file
            - vcf_sample_id: column name in the metafile that indicates the
                sample ID in the VCF file
            - vcf_normal_id: olumn name in the metafile that indicates the
                normal sample ID in the VCF file
            - sex:flagin the metafile that indicates the sample sex
            - guess_baits: The column name in the metafile that indicates
                whether to guess the bait file from the bam files
        baitfile: Potentially targeted genomic regions.
            E.g. all possible exons for the reference genome.
            This is optional when `method` is `wgs`.
        accfile: The accessible genomic regions.
            If not given, use `cnvkit.py access` to generate one.
        access_excludes (list): File(s) with regions to be excluded for
            `cnvkit.py access`.
        guessbaits_guided (flag): Whether to use guided mode for
            guessing baits using `baitfile`, otherwise unguided, using the
            `accfile`.
        guessbaits (flag): Guess the bait file from the bam files,
            either guided or unguided.
            If False, `baitfile` is used. Otherwise, if `baitfile` is given,
            use it (guided), otherwise use `accfile` (unguided).
            The bam files with `metacols.guess_baits` column set to
            `True`, `TRUE`, `true`, `1`, `Yes`, `YES`, or `yes`
            will be used to guess the bait file.
        heatmap_cnr (flag): Whether to generate a heatmap of the
            `.cnr` files (bin-level signals). This is allowed to set to `False`,
            it will take longer to run.
        case: The group name of samples in `metacols.group` to call CNVs for.
            If not specified, use all samples. In such a case, `control` must
            not be specified, as we are using a flat reference.
        control: The group name of samples in `metacols.group` to use as
            reference if not specified, use a flat reference.
        cnvkit: the path to the cnvkit.py executable, defaults to
            `config.exe.cnvkit` from `./.biopipen.toml` or `~/.biopipen.toml`.
        rscript: Path to the Rscript excecutable to use for running R code.
            Requires `DNAcopy` to be installed in R, defaults to
            `config.lang.rscript`
        samtools: Path to samtools, used for guessing bait file.
        convert: Linux `convert` command to convert pdf to png
            So that they can be embedded in the HTML report.
        ncores: Default number of cores to use for all processes with
            `envs.ncores`, defaults to `config.misc.ncores`
        reffa: the reference genome (e.g. hg19.fa).
            Used by `CNVkitAccess`, `CNVkitAutobin` and `CNVkitReference`
        annotate: Use gene models from this file to assign names to the
            target regions. Format: UCSC `refFlat.txt` or `ensFlat.txt` file
            (preferred), or BED, interval list, GFF, or similar.
        short_names (flag): Reduce multi-accession bait labels to
            be short and consistent.
        method (choice): Sequencing protocol, determines whether and how to
            use antitarget bins.
            - hybrid: hybridization capture
            - amplicon: targeted amplicon sequencing
            - wgs: whole genome sequencing
        male_reference (flag): Use or assume a male reference
            (i.e. female samples will have +1 log-CNR of chrX; otherwise
            male samples would have -1 chrX).
            Used by `CNVkitReference`, `CNVkitCall`, `CNVkitHeatmapCns` and
            `CNVkitHeatmapCnr`.
        drop_low_coverage (flag): Drop very-low-coverage bins
            before segmentation to avoid false-positive deletions in
            poor-quality tumor samples. Used by `CNVkitSegment` and `CNVkitCall`
        no_gc (flag): Skip GC correction for
            `cnvkit.py reference/fix`.
        no_edge (flag): Skip edge-effect correction for
            `cnvkit.py reference/fix`.
        no_rmask (flag): Skip RepeatMasker correction for
            `cnvkit.py reference/fix`.
            no_* options are used by `CNVkitReference` and `CNVkitFix`
        min_variant_depth (type=int): Minimum read depth for a SNV to be
            displayed in the b-allele frequency plot.
            Used by `CNVkitSegment` and `CNVkitCall`
        zygosity_freq (type=float): Ignore VCF's genotypes (GT field) and
            instead infer zygosity from allele frequencies.
            Used by `CNVkitSegment` and `CNVkitCall`
    """
    DEFAULTS = Diot(
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
        guessbaits_guided=False,
        male_reference=False,
        drop_low_coverage=False,
        min_variant_depth=20,
        no_gc=False,
        no_edge=False,
        no_rmask=False,
        zygosity_freq=0.25,
        metacols=Diot(
            group="Group",
            purity="Purity",
            snpvcf="SnpVcf",
            bam="Bam",
            vcf_sample_id="VcfSampleId",
            vcf_normal_id="VcfNormalId",
            sex="Sex",
            guess_baits="GuessBaits",
        ),
    )

    @cached_property
    def col(self):
        """Get the column names by self.col.<colname>"""
        return _MetaCol(
            self.opts.get("metacols"),
            self.__class__.DEFAULTS.metacols,
        )

    @ProcGroup.add_proc
    def p_metafile(self):
        """Build MetaFile process"""
        from .misc import File2Proc

        @mark(board_config_hidden=True)
        class MetaFile(File2Proc):
            """Pass by the metafile to the next process.

            When the group argument `metafile` is provided, it will be used
            as the input data, otherwise, this process group should be a
            part of a pipeline, and the metafile will be passed by its
            required processes.
            """
            # Do not require metafile, as we could use the pipeline as part of
            # another pipeline, which can generate a metafile
            # Remember to set the dependency in the pipeline:
            # >>> pipeline.procs.MetaFile.requires = [other_pipeline.procs]
            # where other_pipeline.procs generate the metafile
            if self.opts.metafile:
                input_data = [self.opts.metafile]

        return MetaFile

    @ProcGroup.add_proc
    def p_cnvkit_access(self):
        """Build CNVkitAccess process"""
        if self.opts.get("accfile"):
            from .misc import File2Proc

            @mark(board_config_hidden=True)
            class CNVkitAccess(File2Proc):
                """Pass by the access file to the next process."""
                input_data = [self.opts.accfile]
        else:
            from .cnvkit import CNVkitAccess

            excludes = self.opts.get("excludes", [])
            if not isinstance(excludes, (list, tuple)):
                excludes = [excludes]

            @annotate.format_doc(indent=4)
            class CNVkitAccess(CNVkitAccess):
                """{{Summary}}

                **When group argument `accfile` is provided, the arguments won't
                work. The `accfile` will just be passed by to the next
                process.**

                Envs:
                    cnvkit (pgarg): {{Envs.cnvkit.help | indent: 24}}.
                        Defaults to group argument `cnvkit`.
                    ref (pgarg=reffa): {{Envs.ref.help | indent: 24}}.
                        Defaults group argument `reffa`.
                """
                input_data = [excludes]
                envs = {
                    "cnvkit": self.opts.cnvkit,
                    "ref": self.opts.reffa,
                }

        return CNVkitAccess

    @ProcGroup.add_proc
    def p_cnvkit_guessbaits(self):
        """Build CNVkitGuessBaits process"""
        from .cnvkit import CNVkitGuessBaits

        if (
            not self.opts.guessbaits and
            not is_loading_pipeline("-h", "-h+", "--help", "--help+")
        ):
            return None

        def _guess_baits_bams(ch):
            df = _metadf(_1st(ch))
            if self.col.guess_baits not in df:
                # Use all bams
                return df.loc[:, self.col.bam].tolist()

            # Use only specified
            guess_baits = df[self.col.guess_baits]
            return df.loc[
                (guess_baits == True)  # noqa
                | (guess_baits == "True")
                | (guess_baits == "TRUE")
                | (guess_baits == "true")
                | (guess_baits == "1")
                | (guess_baits == 1)
                | (guess_baits == "yes")
                | (guess_baits == "YES")
                | (guess_baits == "Yes"),
                self.col.bam,
            ].tolist()

        if self.opts.guessbaits_guided:
            if not self.opts.baitfile:
                raise ValueError(
                    "`baitfile` must be specified for guided mode "
                    "to guess baits. See: "
                    "https://cnvkit.readthedocs.io/en/stable/scripts.html"
                )

            @annotate.format_doc(indent=4)
            class CNVkitGuessBaits(CNVkitGuessBaits):
                """{{Summary}}

                Envs:
                    cnvkit (pgarg): {{Envs.cnvkit.help | indent: 24}}.
                        Defaults to group argument `cnvkit`.
                    samtools (pgarg): {{Envs.samtools.help | indent: 24}}.
                        Defaults to group argument `samtools`.
                    ncores (pgarg): {{Envs.ncores.help | indent: 24}}.
                        Defaults to group argument `ncores`.
                    ref (pgarg=reffa): {{Envs.ref.help | indent: 24}}.
                        Defaults to group argument `reffa`.
                    guided (pgarg): {{Envs.guided.help | indent: 24}}.
                        Defaults to group argument `guessbaits_guided`.
                """
                requires = self.p_metafile
                input_data = lambda metafile_ch: tibble(
                    bamfiles=[_guess_baits_bams(metafile_ch)],
                    atfile=self.opts.baitfile,
                )
                envs = {
                    "cnvkit": self.opts.cnvkit,
                    "samtools": self.opts.samtools,
                    "ncores": self.opts.ncores,
                    "ref": self.opts.reffa,
                    "guided": True,
                }
        else:  # unguided
            @annotate.format_doc(indent=4)
            class CNVkitGuessBaits(CNVkitGuessBaits):
                """{{Summary}}

                Envs:
                    cnvkit (pgarg): {{Envs.cnvkit.help | indent: 24}}.
                        Defaults to group argument `cnvkit`.
                    samtools (pgarg): {{Envs.samtools.help | indent: 24}}.
                        Defaults to group argument `samtools`.
                    ncores (pgarg): {{Envs.ncores.help | indent: 24}}.
                        Defaults to group argument `ncores`.
                    ref (pgarg=reffa): {{Envs.ref.help | indent: 24}}.
                        Defaults to group argument `reffa`.
                    guided (pgarg): {{Envs.guided.help | indent: 24}}.
                        Defaults to group argument `guessbaits_guided`.
                """
                requires = self.p_metafile, self.p_cnvkit_access
                input_data = lambda metafile_ch, access_ch: tibble(
                    bamfiles=[_guess_baits_bams(metafile_ch)],
                    accessfile=_1st(access_ch),
                )
                envs = {
                    "cnvkit": self.opts.cnvkit,
                    "samtools": self.opts.samtools,
                    "ncores": self.opts.ncores,
                    "ref": self.opts.reffa,
                    "guided": False,
                }

        return CNVkitGuessBaits

    @ProcGroup.add_proc
    def p_cnvkit_autobin(self):
        """Build CNVkitAutobin process"""
        from .cnvkit import CNVkitAutobin

        @annotate.format_doc(indent=3)
        class CNVkitAutobin(CNVkitAutobin):
            """{{Summary}}

            Envs:
                method (pgarg): {{Envs.method.help | indent: 20}}.
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                ref (pgarg=reffa): {{Envs.ref.help | indent: 20}}.
                    Defaults to group argument `reffa`.
                annotate (pgarg): {{Envs.annotate.help | indent: 20}}.
                    Defaults to group argument `annotate`.
                short_names (pgarg): {{Envs.short_names.help | indent:20}}.
                    Defaults to group argument `short_names`.
            """
            if self.p_cnvkit_guessbaits:
                requires = (
                    self.p_metafile,
                    self.p_cnvkit_access,
                    self.p_cnvkit_guessbaits,
                )
                input_data = lambda ch1, ch2, ch3: tibble(
                    bamfiles=[_metadf(_1st(ch1))[self.col.bam].tolist()],
                    accfile=_1st(ch2),
                    baitfile=(
                        _1st(ch3)
                        if self.opts.guessbaits
                        else self.opts.baitfile
                    ),
                )
            else:
                requires = self.p_metafile, self.p_cnvkit_access
                input_data = lambda ch1, ch2: tibble(
                    bamfiles=[_metadf(_1st(ch1))[self.col.bam].tolist()],
                    accfile=_1st(ch2),
                    baitfile=self.opts.baitfile,
                )
            envs = {
                "cnvkit": self.opts.cnvkit,
                "method": self.opts.method,
                "annotate": self.opts.annotate,
                "short_names": self.opts.short_names,
                "ref": self.opts.reffa,
            }

        return CNVkitAutobin

    def _p_cnvkit_coverage(self, anti: bool):
        """Build CNVkitTargetCoverage and CNVkitAntiTargetCoverage processes"""
        from .cnvkit import CNVkitCoverage

        p = Proc.from_proc(
            CNVkitCoverage,
            name="CNVkitCoverageAnittarget" if anti else "CNVkitCoverageTarget",
            requires=[self.p_metafile, self.p_cnvkit_autobin],
            input_data=lambda ch1, ch2: tibble(
                _metadf(_1st(ch1))[self.col.bam].tolist(),
                target_file=ch2[
                    "antitarget_file" if anti else "target_file"
                ].tolist()[0],
            ),
            envs={
                "cnvkit": self.opts.cnvkit,
                "ncores": self.opts.ncores,
                "ref": self.opts.reffa,
            }
        )
        if anti:
            p.__doc__ = """Build the coverage for the anti-target regions"""
        else:
            p.__doc__ = """Build the coverage for the target regions"""

        p.__doc__ += """

        {{* Summary.long }}

        Envs:
            cnvkit (pgarg): {{Envs.cnvkit.help | indent: 16}}.
                Defaults to group argument `cnvkit`.
            ncores (pgarg): {{Envs.ncores.help | indent: 16}}.
                Defaults to group argument `ncores`.
            ref (pgarg=reffa): {{Envs.ref.help | indent: 16}}.
                Defaults to group argument `reffa`.
        """
        return annotate.format_doc(indent=2)(p)

    @ProcGroup.add_proc
    def p_cnvkit_coverage_target(self):
        """Build CNVkitCoverageTarget process"""
        return self._p_cnvkit_coverage(anti=False)

    @ProcGroup.add_proc
    def p_cnvkit_coverage_antitarget(self):
        """Build CNVkitCoverageAntiTarget process"""
        return self._p_cnvkit_coverage(anti=True)

    @ProcGroup.add_proc
    def p_cnvkit_reference(self):
        """Build CNVkitReference process"""
        from .cnvkit import CNVkitReference

        def _input_data(ch1, ch2, ch3, ch4):
            metadf = _metadf(_1st(ch1))

            if self.opts.control:
                # Use control samples to build reference
                control_masks = metadf[self.col.group] == self.opts.control
                covfiles = [
                    ch2.outfile[control_masks].tolist()
                    + ch3.outfile[control_masks].tolist()
                ]
                target_file = None
                antitarget_file = None
                if self.col.sex in metadf:
                    all_sex = metadf[self.col.sex][control_masks].unique()
                    sample_sex = [None] if len(all_sex) > 1 else all_sex[0]
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

        @annotate.format_doc(indent=3)
        class CNVkitReference(CNVkitReference):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                no_gc (pgarg): {{Envs.no_gc.help | indent: 20}}.
                    Defaults to group argument `no_gc`.
                no_edge (pgarg): {{Envs.no_edge.help | indent: 20}}.
                    Defaults to group argument `no_edge`.
                no_rmask (pgarg): {{Envs.no_rmask.help | indent: 20}}.
                    Defaults to group argument `no_rmask`.
                ref (pgarg=reffa): {{Envs.ref.help | indent: 20}}.
                    Defaults to group argument `reffa`.
                male_reference (pgarg): {{
                    Envs.male_reference.help | indent: 20 }}.
                    Defaults to group argument `male_reference`.
            """
            requires = [
                self.p_metafile,
                self.p_cnvkit_coverage_target,
                self.p_cnvkit_coverage_antitarget,
                self.p_cnvkit_autobin,
            ]
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "no_gc": self.opts.no_gc,
                "no_edge": self.opts.no_edge,
                "no_rmask": self.opts.no_rmask,
                "ref": self.opts.reffa,
                "male_reference": self.opts.male_reference,
            }

        return CNVkitReference

    @ProcGroup.add_proc
    def p_cnvkit_fix(self):
        """Build CNVkitFix process"""
        from .cnvkit import CNVkitFix

        if not self.opts.case and self.opts.control:
            raise ValueError(
                "`case` is not specified, meaning using all samples as cases, "
                "but `control` is specified (we can only use a flat reference "
                "in this case)."
            )

        def _input_data(ch1, ch2, ch3, ch4):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

            return tibble(
                target_file=ch2.outfile[tumor_masks],
                antitarget_file=ch3.outfile[tumor_masks],
                reference=ch4.outfile,
                sample_id=metadf["Sample"][tumor_masks],
            )

        @annotate.format_doc(indent=3)
        class CNVkitFix(CNVkitFix):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                no_gc (pgarg): {{Envs.no_gc.help | indent: 20}}.
                    Defaults to group argument `no_gc`.
                no_edge (pgarg): {{Envs.no_edge.help | indent: 20}}.
                    Defaults to group argument `no_edge`.
                no_rmask (pgarg): {{Envs.no_rmask.help | indent: 20}}.
                    Defaults to group argument `no_rmask`.
            """
            requires = [
                self.p_metafile,
                self.p_cnvkit_coverage_target,
                self.p_cnvkit_coverage_antitarget,
                self.p_cnvkit_reference,
            ]
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "no_gc": self.opts.no_gc,
                "no_edge": self.opts.no_edge,
                "no_rmask": self.opts.no_rmask,
            }

        return CNVkitFix

    @ProcGroup.add_proc
    def p_cnvkit_segment(self):
        """Build CNVkitSegment process"""
        from .cnvkit import CNVkitSegment

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

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

        @annotate.format_doc(indent=3)
        class CNVkitSegment(CNVkitSegment):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                rscript (pgarg): {{Envs.rscript.help | indent: 20}}.
                    Defaults to group argument `rscript`.
                ncores (pgarg): {{Envs.ncores.help | indent: 20}}.
                    Defaults to group argument `ncores`.
                drop_low_coverage (pgarg): {{
                    Envs.drop_low_coverage.help | indent: 20}}.
                    Defaults to group argument `drop_low_coverage`.
                min_variant_depth (pgarg): {{
                    Envs.min_variant_depth.help | indent: 20}}.
                    Defaults to group argument `min_variant_depth`.
                zygosity_freq (pgarg): {{
                    Envs.zygosity_freq.help | indent: 20}}.
                    Defaults to group argument `zygosity_freq`.
            """
            requires = self.p_metafile, self.p_cnvkit_fix
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "rscript": self.opts.rscript,
                "ncores": self.opts.ncores,
                "drop_low_coverage": self.opts.drop_low_coverage,
                "min_variant_depth": self.opts.min_variant_depth,
                "zygosity_freq": self.opts.zygosity_freq,
            }

        return CNVkitSegment

    @ProcGroup.add_proc
    def p_cnvkit_scatter(self):
        """Build CNVkitScatter process"""
        from .cnvkit import CNVkitScatter

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

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

        @annotate.format_doc(indent=3)
        class CNVkitScatter(CNVkitScatter):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                convert (pgarg): {{Envs.convert.help | indent: 20}}.
                    Defaults to group argument `convert`.
                min_variant_depth (pgarg): {{
                    Envs.min_variant_depth.help | indent: 20}}.
                    Defaults to group argument `min_variant_depth`.
            """
            requires = self.p_metafile, self.p_cnvkit_fix, self.p_cnvkit_segment
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "convert": self.opts.convert,
                "min_variant_depth": self.opts.min_variant_depth,
            }

        return CNVkitScatter

    @ProcGroup.add_proc
    def p_cnvkit_diagram(self):
        """Build CNVkitDiagram process"""
        from .cnvkit import CNVkitDiagram

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

            return tibble(
                chrfile=ch2.outfile,
                cnsfile=ch3.outfile,
                sample_sex=(
                    metadf[self.col.sex][tumor_masks]
                    if self.col.sex in metadf
                    else [None]
                ),
            )

        @annotate.format_doc(indent=3)
        class CNVkitDiagram(CNVkitDiagram):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                convert (pgarg): {{Envs.convert.help | indent: 20}}.
                    Defaults to group argument `convert`.
                male_reference (pgarg): {{
                    Envs.male_reference.help | indent: 20}}.
                    Defaults to group argument `male_reference`.
            """
            requires = self.p_metafile, self.p_cnvkit_fix, self.p_cnvkit_segment
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "convert": self.opts.convert,
                "male_reference": self.opts.male_reference,
            }

        return CNVkitDiagram

    @ProcGroup.add_proc
    def p_cnvkit_heatmap_cns(self):
        """Build CNVkitHeatmapCns process"""
        from .cnvkit import CNVkitHeatmap

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

            if self.col.sex in metadf:
                all_sex = metadf[self.col.sex][tumor_masks].unique()
                sample_sex = [None] if len(all_sex) > 1 else all_sex[0]
            else:
                sample_sex = [None]

            return tibble(
                segfiles=[ch2.outfile.tolist()],
                sample_sex=sample_sex,
            )

        @annotate.format_doc(indent=3)
        class CNVkitHeatmapCns(CNVkitHeatmap):
            """Generate heatmaps of segment-level signals of multiple samples

            {{* Summary.long }}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                convert (pgarg): {{Envs.convert.help | indent: 20}}.
                    Defaults to group argument `convert`.
                male_reference (pgarg): {{
                    Envs.male_reference.help | indent: 20}}.
                    Defaults to group argument `male_reference`.
            """
            requires = self.p_metafile, self.p_cnvkit_segment
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "convert": self.opts.convert,
                "male_reference": self.opts.male_reference,
            }

        return CNVkitHeatmapCns

    @ProcGroup.add_proc
    def p_cnvkit_heatmap_cnr(self):
        """Build CNVkitHeatmapCnr process"""
        from .cnvkit import CNVkitHeatmap

        if not self.opts.heatmap_cnr:
            return None

        def _input_data(ch1, ch2):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

            if self.col.sex in metadf:
                all_sex = metadf[self.col.sex][tumor_masks].unique()
                sample_sex = [None] if len(all_sex) > 1 else all_sex[0]
            else:
                sample_sex = [None]

            return tibble(
                segfiles=[ch2.outfile.tolist()],
                sample_sex=sample_sex,
            )

        @annotate.format_doc(indent=3)
        class CNVkitHeatmapCnr(CNVkitHeatmap):
            """Heatmap of bin-level signals of multiple samples

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                convert (pgarg): {{Envs.convert.help | indent: 20}}.
                    Defaults to group argument `convert`.
                male_reference (pgarg): {{
                    Envs.male_reference.help | indent: 20}}.
                    Defaults to group argument `male_reference`.
            """
            requires = self.p_metafile, self.p_cnvkit_fix
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "convert": self.opts.convert,
                "male_reference": self.opts.male_reference,
            }

        return CNVkitHeatmapCnr

    @ProcGroup.add_proc
    def p_cnvkit_call(self):
        """Build CNVkitCall process"""
        from .cnvkit import CNVkitCall

        def _input_data(ch1, ch2, ch3):
            metadf = _metadf(_1st(ch1))
            if not self.opts.case:
                tumor_masks = [True] * len(metadf)
            else:
                tumor_masks = metadf[self.col.group] == self.opts.case

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

        @annotate.format_doc(indent=3)
        class CNVkitCall(CNVkitCall):
            """{{Summary}}

            Envs:
                cnvkit (pgarg): {{Envs.cnvkit.help | indent: 20}}.
                    Defaults to group argument `cnvkit`.
                drop_low_coverage (pgarg): {{
                    Envs.drop_low_coverage.help | indent: 20}}.
                    Defaults to group argument `drop_low_coverage`.
                male_reference (pgarg): {{
                    Envs.male_reference.help | indent: 20}}.
                    Defaults to group argument `male_reference`.
                min_variant_depth (pgarg): {{
                    Envs.min_variant_depth.help | indent: 20}}.
                    Defaults to group argument `min_variant_depth`.
                zygosity_freq (pgarg): {{
                    Envs.zygosity_freq.help | indent: 20}}.
                    Defaults to group argument `zygosity_freq`.
            """
            requires = self.p_metafile, self.p_cnvkit_fix, self.p_cnvkit_segment
            input_data = _input_data
            envs = {
                "cnvkit": self.opts.cnvkit,
                "drop_low_coverage": self.opts.drop_low_coverage,
                "male_reference": self.opts.male_reference,
                "min_variant_depth": self.opts.min_variant_depth,
                "zygosity_freq": self.opts.zygosity_freq,
            }

        return CNVkitCall


if __name__ == "__main__":
    CNVkitPipeline().as_pipen(
        # If we run this procgroup as a whole, we don't want to collapse
        # the processes in the index page of report.
        plugin_opts={"report_no_collapse_pgs": True}
    ).run()
