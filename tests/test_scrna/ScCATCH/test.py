from pathlib import Path
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    CellTypeAnnotation as CellTypeAnnotation_,
    SeuratClusterStats as SeuratClusterStats_,
)
from biopipen.core.testing import get_pipeline


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
        data <- NormalizeData(data)
        data <- FindVariableFeatures(data)
        data <- ScaleData(data)
        data <- RunPCA(data)
        data <- FindNeighbors(data)
        data <- FindClusters(data)

        saveRDS(data, {{out.outfile | quote}})
    """


class CellTypeAnnotation(CellTypeAnnotation_):
    requires = PrepareQuery
    envs = {
        "tool": "sccatch",
        "sccatch_args": {
            "species": "Human",
            "tissue": "Peripheral blood",
        },
    }


class CellTypeAnnotation2(CellTypeAnnotation_):
    requires = PrepareQuery
    envs = {
        "tool": "sccatch",
        "sccatch_args": {
            "marker": str(Path(__file__).parent.joinpath("data", "tcell.sccatch.RDS")),
        },
    }


class SeuratClusterStats(SeuratClusterStats_):
    requires = CellTypeAnnotation
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group-by": "seurat_clusters",
            }
        }
    }


class SeuratClusterStats2(SeuratClusterStats_):
    requires = CellTypeAnnotation2
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group-by": "seurat_clusters",
            }
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareQuery)


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
