"""The splitpipe pipeline

Parse Biosciences pipeline (split-pipe) processes plus a combine process to
merge the sublibraries of a sequencing run and a summary process to summarize
the metrics of multiple (combined) runs.
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

    Args:
        input (list): The list of lists of fastq (R1) files, one list per
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
        fasta (list): The reference sequence file(s) in fasta format, to build
            the genome with `SplitPipeRef`.
        genes (list): The reference gene annotation file(s) in GTF format, to
            build the genome with `SplitPipeRef`.
        genome_name (list): The name(s) of the genome(s), used to name the
            genome directory built by `SplitPipeRef`.
        genome_dir: The directory of an already built genome. If not given,
            the genome directory is built by `SplitPipeRef` from
            `fasta`/`genes`/`genome_name`.
        combine (flag): Whether to combine the sublibraries into a single
            dataset with `SplitPipeCombine` (`--mode comb`).
        summary (flag): Whether to summarize the metrics of the sublibrary
            (or combined) output directories with `SplitPipeSummary`.
    """
    DEFAULTS = Diot(
        input=None,
        fq2=None,
        ids=None,
        fasta=None,
        genes=None,
        genome_name=None,
        genome_dir=None,
        combine=True,
        summary=True,
    )

    def post_init(self):
        """Check the input and the genome"""
        if not is_loading_pipeline("-h", "-h+", "--help", "--help+"):
            if (
                not isinstance(self.opts.input, (list, tuple))
                or len(self.opts.input) == 0
            ):
                raise TypeError(
                    "The input of `SplitPipePipeline` should be a list of lists of "
                    "fastq files, one per sublibrary."
                )
            if not self.opts.genome_dir and not (
                self.opts.fasta and self.opts.genes and self.opts.genome_name
            ):
                raise TypeError(
                    "The `genome_dir` is required unless `fasta`, `genes` and "
                    "`genome_name` are provided, so that the genome can be "
                    "built by `SplitPipeRef`."
                )
            if self.opts.fq2 and (
                not isinstance(self.opts.fq2, (list, tuple))
                or len(self.opts.fq2) != len(self.opts.input)
            ):
                raise TypeError(
                    "The `fq2` of `SplitPipePipeline` should be a list of lists "
                    "of fastq files, one per sublibrary (same length as `input`)."
                )

        if isinstance(self.opts.input, (list, tuple)):
            self.opts.input = [
                [y.strip() for y in x.split(",")]
                if isinstance(x, str)
                else x
                for x in self.opts.input
            ]
        if isinstance(self.opts.fq2, (list, tuple)):
            self.opts.fq2 = [
                [y.strip() for y in x.split(",")]
                if isinstance(x, str)
                else x
                for x in self.opts.fq2
            ]

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_ref(self) -> Type[Proc] | None:
        """Build SplitPipeRef process"""
        if self.opts.genome_dir:
            return None

        from .splitpipe import SplitPipeRef as _SplitPipeRef

        fasta = self.opts.fasta
        if not isinstance(fasta, (list, tuple)):
            fasta = [fasta]
        genes = self.opts.genes
        if not isinstance(genes, (list, tuple)):
            genes = [genes]
        genome_name = self.opts.genome_name

        class SplitPipeRef(_SplitPipeRef):
            input_data = [(fasta, genes)]
            # pass the genome name(s) to name the genome directory and to
            # validate the pairing with the fasta/genes files
            envs = {"genome_name": genome_name}

        return SplitPipeRef

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_count(self) -> Type[Proc]:
        """Build SplitPipeCount process"""
        from .splitpipe import SplitPipeCount as _SplitPipeCount

        # `fq2` may be absent (the R2 names are then derived from R1 by
        # split-pipe), in which case `None` is passed as the input
        inputs_ = self.opts.input or []  # empty when loading the pipeline (--help)
        fq2_ = self.opts.fq2 or [None] * len(inputs_)
        input_ = (
            list(zip(inputs_, fq2_, self.opts.ids))
            if self.opts.ids
            else [(fq1, fq2, None) for fq1, fq2 in zip(inputs_, fq2_)]
        )

        if self.opts.genome_dir:
            genome_dir = self.opts.genome_dir

            class SplitPipeCount(_SplitPipeCount):
                input = "fq1:files, fq2:files, id, genome_dir:dir"
                input_data = [(*fq1fq2id, genome_dir) for fq1fq2id in input_]

            return SplitPipeCount

        class SplitPipeCount(_SplitPipeCount):
            input = "fq1:files, fq2:files, id, genome_dir:dir"
            requires = self.p_splitpipe_ref
            input_data = lambda ch: [
                (*fq1fq2id, ch.iloc[0, 0]) for fq1fq2id in input_
            ]

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
            self.p_splitpipe_combine if self.opts.combine
            else self.p_splitpipe_count
        )

        class SplitPipeSummary(_SplitPipeSummary):
            requires = requires_
            input_data = lambda ch: [list(ch.iloc[:, 0])]

        return SplitPipeSummary
