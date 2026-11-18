"""The splitpipe pipeline

Parse Biosciences pipeline (split-pipe) processes plus a combine process to
merge the sublibraries of a sequencing run and a summary process to summarize
the metrics of multiple (combined) runs.

A sequencing run may also contain a TCR or BCR library for the same
well-pools as the GEX libraries. Pass them as the `tcr`/`bcr` blocks and they
are processed next to the GEX sublibraries (see `SplitPipePipeline`).

The reference genome is not built by this pipeline: build it once with the
`SplitPipeRef` process and pass its output directory as `genome_dir`.
"""

from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING, Type

from argx import Namespace  # type: ignore
from diot import Diot  # type: ignore
from pipen.utils import is_loading_pipeline
from pipen_args.procgroup import ProcGroup

if TYPE_CHECKING:
    from pipen import Proc  # type: ignore


# The immune genomes supported by split-pipe for each analysis
IMMUNE_GENOMES = {
    "tcr": ("human", "mouse"),
    "bcr": ("human", "mouse", "transgenic_mouse"),
}


class SplitPipePipeline(ProcGroup):
    """The splitpipe pipeline

    Run split-pipe for the sublibraries of one sequencing run, combine them
    into a single dataset and summarize the metrics.

    The reference genome must be built beforehand with the `SplitPipeRef`
    process (`split-pipe --mode mkref`), and its output directory passed as
    `genome_dir`.

    A sequencing run may also contain a TCR or BCR library for the same
    well-pools as the GEX sublibraries. Pass it as the `tcr`/`bcr` block: its
    sublibraries are processed by dedicated `SplitPipeTCRCount`/
    `SplitPipeTCRCombine` (resp. `SplitPipeBCRCount`/`SplitPipeBCRCombine`)
    processes. The immune output directories are suffixed with `-tcr`/`-bcr`
    (and `-tcr-comb`/`-bcr-comb` when combined) so they never collide with the
    GEX ones.

    Pairing rule: unless the block sets `parent_dir`, each immune sublibrary
    id must equal the id of a GEX sublibrary, and that GEX output directory
    is used as the immune parent. The immune combination then uses this run's
    combined GEX output as parent, so `combine = true` on the block requires
    `combine = true` on the pipeline. With an explicit `parent_dir` (a GEX
    parent processed by another run), the immune ids are not matched against
    the GEX ids and `combine` does not depend on the pipeline `combine`.

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
        tcr (dict): The TCR sublibrary block, or `None` when the run has no
            TCR library. The block accepts `fq1`, `fq2` and `ids` in the same
            list-of-lists format as the top-level options, plus
            `immune_genome` (`human` by default, or `mouse`), `use_imgt_db`
            (default `False`), `combine` (default `True`) and `parent_dir`
            (the GEX parent directory, when it was produced by a different
            run). Unless `parent_dir` is given, each immune id must equal a
            GEX sublibrary id.
        bcr (dict): The BCR sublibrary block, in the same format as `tcr`,
            except that `immune_genome` can also be `transgenic_mouse`.
    """

    DEFAULTS = Diot(
        fq1=None,
        fq2=None,
        ids=None,
        genome_dir=None,
        combine=True,
        summary=True,
        tcr=None,
        bcr=None,
    )

    def post_init(self):
        """Check the input and the genome"""
        loading = is_loading_pipeline("-h", "-h+", "--help", "--help+")
        if not loading:
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

        for kind in ("tcr", "bcr"):
            if self.opts.get(kind) is not None:
                self._check_immune_block(kind, loading)

    def _check_immune_block(self, kind: str, loading: bool) -> None:
        """Normalize and validate a `tcr`/`bcr` block"""
        block = self.opts[kind]
        if isinstance(block, Namespace):
            block = vars(block)
        if not isinstance(block, dict):
            raise TypeError(
                f"`{kind}` of `SplitPipePipeline` should be a dict with the "
                f"fastq files, ids and immune options of the {kind.upper()} "
                f"sublibraries."
            )
        block = Diot(block)
        self.opts[kind] = block

        if isinstance(block.get("fq1"), (list, tuple)):
            block["fq1"] = [
                [y.strip() for y in x.split(",")] if isinstance(x, str) else x
                for x in block["fq1"]
            ]
        if isinstance(block.get("fq2"), (list, tuple)):
            block["fq2"] = [
                [y.strip() for y in x.split(",")] if isinstance(x, str) else x
                for x in block["fq2"]
            ]

        block.setdefault("fq2", None)
        block.setdefault("ids", None)
        block.setdefault("parent_dir", None)
        block.setdefault("combine", True)
        block.setdefault("use_imgt_db", False)
        block.setdefault("immune_genome", "human")

        if loading:
            return

        if (
            not isinstance(block.fq1, (list, tuple))
            or not block.fq1
            or not all(isinstance(fq1, (list, tuple)) for fq1 in block.fq1)
        ):
            raise TypeError(
                f"`{kind}.fq1` of `SplitPipePipeline` should be a non-empty "
                f"list of lists of fastq files, one per {kind.upper()} "
                f"sublibrary."
            )
        if block.fq2 is not None and (
            not isinstance(block.fq2, (list, tuple))
            or len(block.fq2) != len(block.fq1)
        ):
            raise TypeError(
                f"`{kind}.fq2` of `SplitPipePipeline` should be a list of "
                f"lists of fastq files, one per {kind.upper()} sublibrary "
                f"(same length as `{kind}.fq1`)."
            )
        if not isinstance(block.ids, (list, tuple)) or len(block.ids) != len(
            block.fq1
        ):
            raise TypeError(
                f"`{kind}.ids` of `SplitPipePipeline` should be a list of "
                f"ids, one per {kind.upper()} sublibrary (same length as "
                f"`{kind}.fq1`)."
            )
        if block.immune_genome not in IMMUNE_GENOMES[kind]:
            raise ValueError(
                f"`{kind}.immune_genome` should be one of "
                f"{IMMUNE_GENOMES[kind]}; got {block.immune_genome!r}."
            )

        if block.parent_dir is not None:
            return

        if not self.opts.ids:
            raise TypeError(
                f"`ids` of `SplitPipePipeline` is required to pair the "
                f"{kind.upper()} sublibraries with the GEX ones. Set it, or "
                f"set `{kind}.parent_dir` when the GEX parent was processed "
                f"by another run."
            )
        unmatched = [id_ for id_ in block.ids if id_ not in self.opts.ids]
        if unmatched:
            raise ValueError(
                f"The {kind.upper()} sublibrary ids {unmatched} are not GEX "
                f"sublibrary ids {list(self.opts.ids)}. The {kind.upper()} id "
                f"must equal its GEX sublibrary id, or set `{kind}.parent_dir` "
                f"when the GEX parent was processed by another run."
            )
        if block.combine and not self.opts.combine:
            raise ValueError(
                f"`{kind}.combine = true` requires `combine = true` on "
                f"`SplitPipePipeline`, so the {kind.upper()} sublibraries can "
                f"be combined against the combined GEX output. Set "
                f"`combine = true`, or `{kind}.combine = false`, or set "
                f"`{kind}.parent_dir` when the GEX parent was processed by "
                f"another run."
            )

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
            input_data = input_
            envs = {"genome_dir": self.opts.genome_dir}

        return SplitPipeCount

    def _immune_input_data(self, kind: str):
        """Build the input_data of the immune count process"""
        block = self.opts[kind]
        fq2_ = block.fq2 or [None] * len(block.fq1)

        if block.parent_dir is not None:
            rows = [
                (fq1, fq2, id_, self.opts.genome_dir, block.parent_dir)
                for fq1, fq2, id_ in zip(block.fq1, fq2_, block.ids)
            ]

            def input_data(ch):
                return rows

            return input_data

        def input_data(ch):
            gex_dirs = {Path(str(outdir)).name: outdir for outdir in ch.iloc[:, 0]}
            unmatched = [id_ for id_ in block.ids if id_ not in gex_dirs]
            if unmatched:
                raise ValueError(
                    f"The {kind.upper()} sublibrary ids {unmatched} do not "
                    f"match any GEX sublibrary output directory. The "
                    f"{kind.upper()} id must equal its GEX sublibrary id, or "
                    f"set `{kind}.parent_dir`."
                )
            return [
                (fq1, fq2, id_, self.opts.genome_dir, gex_dirs[id_])
                for fq1, fq2, id_ in zip(block.fq1, fq2_, block.ids)
            ]

        return input_data

    def _build_immune_count(self, kind: str) -> Type[Proc] | None:
        """Build the SplitPipeCount process of a `tcr`/`bcr` block"""
        block = self.opts.get(kind)
        if block is None:
            return None

        from .splitpipe import SplitPipeCount as _SplitPipeCount

        input_data_ = self._immune_input_data(kind)
        envs_ = {
            "genome_dir": self.opts.genome_dir,
            f"{kind}_analysis": True,
            "immune_genome": block.immune_genome,
            "use_imgt_db": block.use_imgt_db,
        }
        if block.parent_dir is not None:
            envs_["parent_dir"] = block.parent_dir

        class SplitPipeCount(_SplitPipeCount):
            name = f"SplitPipe{kind.upper()}Count"
            input = "fq1:files, fq2:files, id, genome_dir:dir, parent_dir:dir"
            input_data = input_data_
            output = f"outdir:dir:{{{{in.id}}}}-{kind}"
            envs = envs_
            requires = self.p_splitpipe_count

        return SplitPipeCount

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_tcr_count(self) -> Type[Proc] | None:
        """Build the SplitPipeCount process for the TCR sublibraries"""
        return self._build_immune_count("tcr")

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_bcr_count(self) -> Type[Proc] | None:
        """Build the SplitPipeCount process for the BCR sublibraries"""
        return self._build_immune_count("bcr")

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

    def _build_immune_combine(self, kind: str) -> Type[Proc] | None:
        """Build the SplitPipeCombine process of a `tcr`/`bcr` block"""
        block = self.opts.get(kind)
        if block is None or not block.combine:
            return None

        from .splitpipe import SplitPipeCombine as _SplitPipeCombine

        immune_count = getattr(self, f"p_splitpipe_{kind}_count")
        envs_ = {
            "immune_genome": block.immune_genome,
            "use_imgt_db": block.use_imgt_db,
        }

        if block.parent_dir is not None:
            requires_ = immune_count
            envs_["parent_dir"] = block.parent_dir

            def input_data_(ch):
                return [(list(ch.iloc[:, 0]), None, str(block.parent_dir))]

        else:
            requires_ = [immune_count, self.p_splitpipe_combine]

            def input_data_(ch_immune, ch_gex):
                if (
                    ch_immune.shape[0] < 1
                    or ch_immune.shape[1] < 1
                    or ch_gex.shape[0] < 1
                    or ch_gex.shape[1] < 1
                ):
                    raise ValueError(
                        f"Cannot combine the {kind.upper()} sublibraries: "
                        f"expected the immune and GEX combine channels, got "
                        f"shapes {ch_immune.shape} and {ch_gex.shape}."
                    )
                gex_comb_dir = ch_gex.iloc[0, 0]
                if not gex_comb_dir:
                    raise ValueError(
                        f"The combined GEX output directory is missing; "
                        f"cannot set `parent_dir` for the {kind.upper()} "
                        f"combination."
                    )
                return [(list(ch_immune.iloc[:, 0]), None, str(gex_comb_dir))]

        class SplitPipeCombine(_SplitPipeCombine):
            name = f"SplitPipe{kind.upper()}Combine"
            input = "sublibs:dirs, id, parent_dir:dir"
            requires = requires_
            envs = envs_
            input_data = input_data_

        return SplitPipeCombine

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_tcr_combine(self) -> Type[Proc] | None:
        """Build the SplitPipeCombine process for the TCR sublibraries"""
        return self._build_immune_combine("tcr")

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_bcr_combine(self) -> Type[Proc] | None:
        """Build the SplitPipeCombine process for the BCR sublibraries"""
        return self._build_immune_combine("bcr")

    @ProcGroup.add_proc  # type: ignore
    def p_splitpipe_summary(self) -> Type[Proc] | None:
        """Build SplitPipeSummary process"""
        if not self.opts.summary:
            return None

        from .splitpipe import SplitPipeSummary as _SplitPipeSummary

        requires_ = [
            self.p_splitpipe_combine if self.opts.combine else self.p_splitpipe_count
        ]
        for kind in ("tcr", "bcr"):
            block = self.opts.get(kind)
            if block is None:
                continue
            requires_.append(
                getattr(
                    self,
                    f"p_splitpipe_{kind}_{'combine' if block.combine else 'count'}",
                )
            )

        if len(requires_) == 1:
            requires_ = requires_[0]

            def input_data_(ch):
                return [list(ch.iloc[:, 0])]

        else:

            def input_data_(*chs):
                return [[x for ch in chs for x in ch.iloc[:, 0]]]

        class SplitPipeSummary(_SplitPipeSummary):
            requires = requires_
            input_data = input_data_

        return SplitPipeSummary
