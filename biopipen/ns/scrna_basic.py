"""Basic analysis for single cell RNA-seq data

- QC
- Clustering
- Marker genes
- Enrichment analysis
"""
from __future__ import annotations
from pathlib import Path
from typing import Type

from pipen.utils import mark
from pipen_annotate import annotate
from pipen_args import ProcGroup
from pipen_board import from_pipen_board

from ..core.proc import Proc


class ScrnaBasic(ProcGroup):
    """Basic analysis for single cell RNA-seq data

    Including QC, clustering, marker genes, and enrichment analysis.

    Args:
        infile: The input file. Either a tab-delimited file containing
            the information of metadata and paths to results of cellranger
            or a seurat object has been saved as RDS file (with extension
            `.rds` or `.RDS`), which QC is assumed to be done.
            As for the tab-delimited file, it should have two columns:
            `Sample` and `RNAData`. `Sample` should be the first column with
            unique identifiers for the samples and `RNAData` indicates where the
            barcodes, genes, expression matrices are.
        is_seurat (flag): Whether the input file is a seurat object
            in RDS format.
            If this process group runs independently, this argument should
            not be set. It will be recognized automatically by the extension
            of `infile`. However, if this process group is run as a part of
            a pipeline, this argument should be set manually since `infile`
            should not be set in this case. It will be passed by other processes
        clustering (choice;required): Which clustering method to use.
            - supervised: Mapping the cells to given reference.
                Using Seurat Reference Mapping procedure.
                See: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
            - unsupervised: Clustering the cells without reference.
                Using Seurat FindClusters procedure.
            - both: Both supervised and unsupervised clustering.
                Performing both of the above procedures. The unsupervised
                clustering will be added as `seurat_clusters_unsupervised`
                to the metadata.
        ref: The reference file for supervised clustering. It should be an
            RDS file (with extension `.rds` or `.RDS`) containing a seurat
            object, or a h5 file (with extension `.h5` or `.h5seurat`) that
            can be loaded by `Seurat::LoadH5Seurat()`.
    """  # noqa: E501

    DEFAULTS = {
        "infile": None,
        "is_seurat": False,
        "clustering": None,
        "ref": None,
    }

    def post_init(self) -> None:
        if self.opts.infile:
            suffix = Path(self.opts.infile).suffix
            self.opts.is_seurat = suffix in (".rds", ".RDS")

    @ProcGroup.add_proc
    def p_input(self) -> Type[Proc]:
        """Build the input for the process group"""
        from .misc import File2Proc

        @mark(board_config_hidden=True)
        class ScrnaBasicInput(File2Proc):
            """Input file for scrna_basic process group

            To specify the input file, use the `infile` argument of the
            process group.
            """

            if self.opts.infile:
                input_data = [self.opts.infile]

        return ScrnaBasicInput

    @ProcGroup.add_proc
    def p_prepare(self) -> Type[Proc]:
        """Prepare the input data into a Seurat object and do QC"""
        if self.opts.is_seurat:
            return self.p_input

        from .scrna import SeuratPreparing

        class ScrnaBasicPrepareAndQC(SeuratPreparing):
            requires = self.p_input

        return ScrnaBasicPrepareAndQC

    @ProcGroup.add_proc
    def p_supervised(self) -> Type[Proc]:
        if (
            self.opts.clustering == "unsupervised"
            and not from_pipen_board()
        ):
            return None

        from .scrna import SeuratMap2Ref

        @annotate.format_doc(indent=3)
        class ScrnaBasicSupervised(SeuratMap2Ref):
            """{{Summary}}

            **Only available when the group argument `clustering` is set to
            `supervised` or `both`.**

            Envs:
                ref (pgarg): {{Envs.ref.help | indent(20)}}.
                    Defaults to the `ref` argument of the process group.
            """
            requires = self.p_prepare
            envs = {
                "ref": self.opts.ref,
            }

        return ScrnaBasicSupervised

    @ProcGroup.add_proc
    def p_supervised_stats(self) -> Type[Proc]:
        if not self.p_supervised and not from_pipen_board():
            return None

        from .scrna import SeuratClusterStats

        @annotate.format_doc(indent=3)
        class ScrnaBasicSupervisedStats(SeuratClusterStats):
            """{{Summary}}

            **Only available when the group argument `clustering` is set to
            `supervised` or `both`.**
            """
            requires = self.p_supervised

        return ScrnaBasicSupervisedStats

    @ProcGroup.add_proc
    def p_unsupervised(self) -> Type[Proc]:
        if (
            self.opts.clustering == "supervised"
            and not from_pipen_board()
        ):
            return None

        from .scrna import SeuratClustering

        class ScrnaBasicUnsupervised(SeuratClustering):
            requires = self.p_prepare

        return ScrnaBasicUnsupervised

    @ProcGroup.add_proc
    def p_unsupervised_anno(self) -> Type[Proc]:
        if not self.p_unsupervised and not from_pipen_board():
            return None

        from .scrna import CellTypeAnnotation

        class ScrnaBasicAnnotation(CellTypeAnnotation):
            requires = self.p_unsupervised

        return ScrnaBasicAnnotation

    @ProcGroup.add_proc
    def p_unsupervised_stats(self) -> Type[Proc]:
        if not self.p_unsupervised_anno and not from_pipen_board():
            return None

        from .scrna import SeuratClusterStats

        class ScrnaBasicUnsupervisedStats(SeuratClusterStats):
            requires = self.p_unsupervised_anno

        return ScrnaBasicUnsupervisedStats

    @ProcGroup.add_proc
    def p_merge(self) -> Type[Proc]:
        if self.opts.clustering == "supervised" and not from_pipen_board():
            return self.p_supervised

        if self.opts.clustering == "unsupervised" and not from_pipen_board():
            return self.p_unsupervised_anno

        @mark(board_config_hidden=True)
        class ScrnaBasicMerge(Proc):
            """Merge the supervised and unsupervised clustering results

            Add unsupervised clustering as metadata to the seurat object
            with supervised clustering.

            The unsupervised clustering results are stored in the metadata
            `seurat_clusters_unsupervised`.

            **Only available when the group argument `clustering` is set to
            `both`.**
            """
            requires = [self.p_supervised, self.p_unsupervised_anno]
            lang = self.p_supervised.lang
            input = "sobjfile:file, uobjfile:file"
            output = "outfile:file:{{in.sobjfile | stem}}.rds"
            script = """
                library(Seurat)
                sobj <- readRDS({{in.sobjfile | quote}})
                uobj <- readRDS({{in.uobjfile | quote}})
                umeta <- as.list(uobj$seurat_clusters)
                names(umeta) <- rownames(uobj)
                sobj <- AddMetaData(
                    sobj,
                    metadata=umeta,
                    col.name="seurat_clusters_unsupervised"
                )
                saveRDS(sobj, {{out.outfile | quote}})
            """

        return ScrnaBasicMerge

    @ProcGroup.add_proc
    def p_findmarkers(self) -> Type[Proc]:
        from .scrna import MarkersFinder

        @annotate.format_doc(indent=3)
        class ScrnaBasicMarkers(MarkersFinder):
            """{{Summary}}

            If the group argument `clustering` is set to `"both"`,
            you can set `group-by` to `"seurat_clusters_unsupervised"` in
            a different case to find the markers for the unsupervised clusters.
            """
            requires = self.p_merge

        return ScrnaBasicMarkers

    @ProcGroup.add_proc
    def p_scgsea(self) -> Type[Proc]:
        from .scrna import ScFGSEA

        class ScrnaBasicScGSEA(ScFGSEA):
            requires = self.p_merge

        return ScrnaBasicScGSEA


if __name__ == "__main__":
    ScrnaBasic().as_pipen().run()
