"""The splitpipe pipeline

Parse Biosciences pipeline (split-pipe) processes plus a combine process to
merge the sublibraries of a sequencing run and a summary process to summarize
the metrics of multiple (combined) runs.

The reference genome is not built by this pipeline: build it once with the
`SplitPipeRef` process and pass its output directory as `genome_dir`.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Type

from diot import Diot  # type: ignore
from pipen.utils import is_loading_pipeline
from pipen_args.procgroup import ProcGroup

if TYPE_CHECKING:
    from pipen import Proc  # type: ignore


class SplitPipePipeline(ProcGroup):
    """The splitpipe pipeline

    Run split-pipe for the sublibraries of one sequencing run, combine them
    into a single dataset and summarize the metrics.

    The reference genome must be built beforehand with the `SplitPipeRef`
    process (`split-pipe --mode mkref`), and its output directory passed as
    `genome_dir`.

    Args:
        fq1 (list): The list of lists of fastq (R1) files, one list per
            sublibrary, or the list of comma-separated strings of fastq files.
            Fastq files of multiple sequencing lanes of the same sublibrary
            belong in the same list.
        fq2 (list): The list of lists of fastq (R2) files, one list per
            sublibrary, or the list of comma-separated strings of fastq files.
            Fastq files of multiple sequencing lanes of the same sublibrary
            belong in the same list. If not given, the R2 file names are
            derived from the R1 ones by split-pipe, which requires the R1 and
            R2 files to differ only by `R1`/`R2` in their names.
        ids (list): The list of ids for the sublibraries.
        genome_dir: The directory of a genome built by `SplitPipeRef`
            (`split-pipe --mode mkref`). Required: run the `SplitPipeRef`
            process once beforehand (see
            `docs/pipelines/splitpipe_pipeline.md`).
        combine (flag): Whether to combine the sublibraries into a single
            dataset with `SplitPipeCombine` (`--mode comb`).
        summary (flag): Whether to summarize the metrics of the sublibrary
            (or combined) output directories with `SplitPipeSummary`.
    """

    DEFAULTS = Diot(
        fq1=None,
        fq2=None,
        ids=None,
        genome_dir=None,
        combine=True,
        summary=True,
    )

    def post_init(self):
        """Check the input and the genome"""
        if not is_loading_pipeline("-h", "-h+", "--help", "--help+"):
            if not isinstance(self.opts.fq1, (list, tuple)) or len(self.opts.fq1) == 0:
                raise TypeError(
                    "The input of `SplitPipePipeline` should be a list of lists of "
                    "fastq files, one per sublibrary."
                )
            if not self.opts.genome_dir:
                raise TypeError(
                    "The `genome_dir` of `SplitPipePipeline` is required: build "
                    "the genome first with the `SplitPipeRef` process (see "
                    "`docs/pipelines/splitpipe_pipeline.md`)."
                )
            if self.opts.fq2 and (
                not isinstance(self.opts.fq2, (list, tuple))
                or len(self.opts.fq2) != len(self.opts.fq1)
            ):
                raise TypeError(
                    "The `fq2` of `SplitPipePipeline` should be a list of lists "
                    "of fastq files, one per sublibrary (same length as `fq1`)."
                )

        if isinstance(self.opts.fq1, (list, tuple)):
            self.opts.fq1 = [
                [y.strip() for y in x.split(",")] if isinstance(x, str) else x
                for x in self.opts.fq1
            ]
        if isinstance(self.opts.fq2, (list, tuple)):
            self.opts.fq2 = [
                [y.strip() for y in x.split(",")] if isinstance(x, str) else x
                for x in self.opts.fq2
            ]

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_count(self) -> Type[Proc]:
        """Build SplitPipeCount process"""
        from .splitpipe import SplitPipeCount as _SplitPipeCount

        # `fq2` may be absent (the R2 names are then derived from R1 by
        # split-pipe), in which case `None` is passed as the input
        inputs_ = self.opts.fq1 or []  # empty when loading the pipeline (--help)
        fq2_ = self.opts.fq2 or [None] * len(inputs_)
        input_ = (
            list(zip(inputs_, fq2_, self.opts.ids))
            if self.opts.ids
            else [(fq1, fq2, None) for fq1, fq2 in zip(inputs_, fq2_)]
        )

        class SplitPipeCount(_SplitPipeCount):
            input = "fq1:files, fq2:files, id, genome_dir:dir"
            input_data = [(*fq1fq2id, self.opts.genome_dir) for fq1fq2id in input_]

        return SplitPipeCount

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_combine(self) -> Type[Proc] | None:
        """Build SplitPipeCombine process"""
        if not self.opts.combine:
            return None

        from .splitpipe import SplitPipeCombine as _SplitPipeCombine

        class SplitPipeCombine(_SplitPipeCombine):
            requires = self.p_splitpipe_count
            input_data = lambda ch: [(list(ch.iloc[:, 0]), None)]

        return SplitPipeCombine

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_summary(self) -> Type[Proc] | None:
        """Build SplitPipeSummary process"""
        if not self.opts.summary:
            return None

        from .splitpipe import SplitPipeSummary as _SplitPipeSummary

        requires_ = (
            self.p_splitpipe_combine if self.opts.combine else self.p_splitpipe_count
        )

        class SplitPipeSummary(_SplitPipeSummary):
            requires = requires_
            input_data = lambda ch: [list(ch.iloc[:, 0])]

        return SplitPipeSummary
