"""The cellranger pipelines

Primarily cellranger process plus summary for summarizing the metrics for
multiple samples.
"""
from __future__ import annotations
from typing import TYPE_CHECKING

from diot import Diot
from pipen.utils import is_loading_pipeline
from pipen_args.procgroup import ProcGroup

if TYPE_CHECKING:
    from pipen import Proc


class CellRangerCountPipeline(ProcGroup):
    """The cellranger count pipeline

    Run cellranger count for multiple samples and summarize the metrics.

    Args:
        input (list): The list of lists of fastq files.
            or the list of comma-separated string of fastq files.
        ids (list): The list of ids for the samples.
    """
    DEFAULTS = Diot(input=None, ids=None)

    def post_init(self):
        """Check if the input is a list of fastq files"""
        if not is_loading_pipeline("-h", "-h+", "--help", "--help+") and (
            not isinstance(self.opts.input, (list, tuple))
            or len(self.opts.input) == 0
        ):
            raise TypeError(
                "The input of `CellRangerCountPipeline` should be a list of lists of "
                "fastq files."
            )

        if isinstance(self.opts.input, (list, tuple)):
            self.opts.input = [
                [y.strip() for y in x.split(",")]
                if isinstance(x, str)
                else x
                for x in self.opts.input
            ]

    @ProcGroup.add_proc
    def p_cellranger_count(self) -> Proc:
        """Build CellRangerCount process"""
        from .cellranger import CellRangerCount as _CellRangerCount

        class CellRangerCount(_CellRangerCount):
            if self.opts.ids:
                input_data = list(zip(self.opts.input, self.opts.ids))
            else:
                input_data = self.opts.input

        return CellRangerCount

    @ProcGroup.add_proc
    def p_cellranger_count_summary(self) -> Proc:
        """Build CellRangerCountSummary process"""
        from .cellranger import CellRangerSummary

        class CellRangerCountSummary(CellRangerSummary):
            requires = self.p_cellranger_count
            input_data = lambda ch: [list(ch.iloc[:, 0])]

        return CellRangerCountSummary


class CellRangerVdjPipeline(ProcGroup):
    """The cellranger vdj pipeline

    Run cellranger vdj for multiple samples and summarize the metrics.

    Args:
        input (list): The list of lists of fastq files.
            or the list of comma-separated string of fastq files.
        ids (list): The list of ids for the samples.
    """
    DEFAULTS = Diot(input=None, ids=None)

    def post_init(self):
        """Check if the input is a list of fastq files"""
        if not is_loading_pipeline("-h", "-h+", "--help", "--help+") and (
            not isinstance(self.opts.input, (list, tuple))
            or len(self.opts.input) == 0
        ):
            raise TypeError(
                "The input of `CellRangerVdjPipeline` should be a list of lists of "
                "fastq files."
            )

        if isinstance(self.opts.input, (list, tuple)):
            self.opts.input = [
                [y.strip() for y in x.split(",")]
                if isinstance(x, str)
                else x
                for x in self.opts.input
            ]

    @ProcGroup.add_proc
    def p_cellranger_vdj(self) -> Proc:
        """Build CellRangerVdj process"""
        from .cellranger import CellRangerVdj as _CellRangerVdj

        class CellRangerVdj(_CellRangerVdj):
            if self.opts.ids:
                input_data = list(zip(self.opts.input, self.opts.ids))
            else:
                input_data = self.opts.input

        return CellRangerVdj

    @ProcGroup.add_proc
    def p_cellranger_vdj_summary(self) -> Proc:
        """Build CellRangerVdjSummary process"""
        from .cellranger import CellRangerSummary

        class CellRangerVdjSummary(CellRangerSummary):
            requires = self.p_cellranger_vdj
            input_data = lambda ch: [list(ch.iloc[:, 0])]

        return CellRangerVdjSummary
