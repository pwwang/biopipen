from pathlib import Path

from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratMap2Ref as SeuratMap2Ref_,
    SeuratClusterStats as SeuratClusterStats_,
    CellTypeAnnotation as CellTypeAnnotation_,
)
from biopipen.core.testing import get_pipeline, _get_test_dirs

_NAME, _WORKDIR, _ = _get_test_dirs(__file__, False)

PBMC3K_REF = str(
    Path(_WORKDIR) / _NAME / "PrepareQuery" / "0" / "output"
    / "pbmc3k.RDS"
)


class PrepareQuery(Proc):
    """Prepare the query data"""
    input = "name"
    input_data = ["pbmc3k"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        library(Seurat)
        library(SeuratData)
        name <- {{in.name | r}}
        InstallData(name)
        data <- LoadData(name)
        data <- UpdateSeuratObject(data)
        data$Sample <- paste0("S", sample(1:2, nrow(data), replace = TRUE))
        data <- NormalizeData(data)
        data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        data <- ScaleData(data, features = rownames(data))
        data <- RunPCA(data, npcs = 30)
        data <- RunUMAP(
            data,
            dims = 1:30,
            n.neighbors = 30,
            min.dist = 0.3,
            return.model = TRUE
        )
        saveRDS(data, {{out.outfile | quote}})
    """


class SeuratMap2Ref(SeuratMap2Ref_):
    requires = PrepareQuery
    order = -1
    envs = {
        "ncores": 2,
        "use": "seurat_annotations",
        "ref": PBMC3K_REF,
        "MapQuery": {"reference.reduction": "pca", "reduction.model": "umap"},
    }


class SeuratMap2RefSplitby(SeuratMap2Ref_):
    requires = PrepareQuery
    envs = {
        "ncores": 2,
        "split_by": "Sample",
        "use": "seurat_annotations",
        "ref": PBMC3K_REF,
        "MapQuery": {"reference.reduction": "pca", "reduction.model": "umap"},
    }


class CellTypeAnnotationMapQuery(CellTypeAnnotation_):
    """The mapquery tool of CellTypeAnnotation on the same query and reference
    as SeuratMap2Ref, cell-level: the reference's seurat_annotations labels are
    transferred to the query cells."""

    requires = PrepareQuery
    envs = {
        "tool": "mapquery",
        "mapquery": {
            "db": PBMC3K_REF,
            "use": "seurat_annotations",
            "MapQuery": {"reference.reduction": "pca", "reduction.model": "umap"},
        },
    }


class CellTypeAnnotationScmap(CellTypeAnnotation_):
    """The scmap tool of CellTypeAnnotation on the same query and reference as
    SeuratMap2Ref, cell-level: the reference's seurat_annotations labels are
    projected onto the query cells."""

    requires = PrepareQuery
    envs = {
        "tool": "scmap",
        "scmap": {"db": PBMC3K_REF, "cluster_col": "seurat_annotations"},
    }


class CellTypeAnnotationCheetah(CellTypeAnnotation_):
    """The cheetah tool of CellTypeAnnotation on the same query and reference as
    SeuratMap2Ref, cell-level: the reference's seurat_annotations labels are
    transferred to the query cells."""

    requires = PrepareQuery
    envs = {
        "tool": "cheetah",
        "cheetah": {"db": PBMC3K_REF, "ref_ct": "seurat_annotations"},
    }


class SeuratClusterStats(SeuratClusterStats_):
    requires = SeuratMap2Ref
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group_by": "seurat_clusters",
            }
        }
    }


class SeuratClusterStats2(SeuratClusterStats_):
    requires = SeuratMap2RefSplitby
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group_by": "seurat_clusters",
            }
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareQuery)
    # return get_pipeline(__file__, enable_report=True).set_starts(PrepareQuery)


def testing(pipen):
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc3k.cluster_stats",
            "dimplots",
            "Dimensional-reduction-plot.dim.png",
        )
    )
    assert outfile.is_file(), str(outfile)

    # CellTypeAnnotation with the reference-based tools: one transferred label
    # per cell
    for name, tool in (
        ("CellTypeAnnotationMapQuery", "mapquery"),
        ("CellTypeAnnotationScmap", "scmap"),
        ("CellTypeAnnotationCheetah", "cheetah"),
    ):
        proc = [proc for proc in pipen.procs if proc.name == name][0]
        outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
        cell_tsv = outfile.with_name(outfile.name + ".cell2celltype.tsv")
        lines = cell_tsv.read_text().splitlines()
        assert lines[0] == "Cell\tDEFAULT"
        labels = [line.split("\t")[1] for line in lines[1:]]
        assert len(labels) > 1 and len(set(labels)) > 1 and all(labels)
        assert (
            f"Renamed annotation column '{tool}_celltype' to 'CellType'"
            in proc.workdir.joinpath("0", "job.stdout").read_text()
        )


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
