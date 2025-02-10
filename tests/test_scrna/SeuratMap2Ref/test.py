from pathlib import Path

from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratMap2Ref as SeuratMap2Ref_,
    SeuratClusterStats,
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
        data$Sample <- paste0("S", sample(1:2, nrow(data), replace = TRUE))
        saveRDS(data, {{out.outfile | quote}})
    """


class SeuratMap2Ref(SeuratMap2Ref_):
    requires = PrepareQuery
    envs = {
        "ncores": 2,
        "use": "celltype.l2",
        "ref": str(
            Path(__file__).parent.parent.parent
            / "data"
            / "reference"
            / "pbmc_multimodal_2023.rds"
        ),
    }


class SeuratMap2Ref2(SeuratMap2Ref_):
    requires = PrepareQuery
    envs = {
        "ncores": 2,
        "split_by": "Sample",
        "use": "celltype.l2",
        "ref": str(
            Path(__file__).parent.parent.parent
            / "data"
            / "reference"
            / "pbmc_multimodal_2023.rds"
        ),
    }


class SeuratClusterStats(SeuratClusterStats):
    requires = SeuratMap2Ref
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group-by": "seurat_clusters",
            }
        }
    }


class SeuratClusterStats2(SeuratClusterStats):
    requires = SeuratMap2Ref2
    envs = {
        "stats": {
            "Number of cells in each cluster by Sample": {
                "group-by": "seurat_clusters",
            }
        }
    }


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_starts(PrepareQuery)
    )


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


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
