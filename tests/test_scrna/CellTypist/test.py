from pathlib import Path

from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    CellTypeAnnotation,
    AnnData2Seurat as AnnData2Seurat_,
    Seurat2AnnData as Seurat2AnnData_,
)
from biopipen.ns.web import Download
from biopipen.core.testing import get_pipeline

MODEL = Path(__file__).parent / "data" / "Immune_All_High.pkl"


class DownloadDemoData(Download):
    input_data = [
        "https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_2000_cells.h5ad"
    ]
    envs = {"tool": "aria2c"}


class AnnData2Seurat(AnnData2Seurat_):
    requires = DownloadDemoData
    envs = {"ident": "cell_type"}


class ModifyCellType(Proc):
    requires = AnnData2Seurat
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.qs"
    lang = config.lang.rscript
    script = """
    library(biopipen.utils)
    library(Seurat)
    sobjfile <- {{in.infile | r}}
    outfile <- {{out.outfile | r}}
    sobj <- read_obj(sobjfile)
    sobj$cell_type <- factor(paste0("c", as.integer(as.factor(sobj$cell_type)) + 1))
    Idents(sobj) <- "cell_type"
    write_obj(sobj, outfile)
    """


class CellTypeAnnotationAnnData(CellTypeAnnotation):
    requires = DownloadDemoData
    envs = {
        "tool": "celltypist",
        "celltypist_args": {"model": MODEL},
    }


class CellTypeAnnotationAnnDataNoOverClustering(CellTypeAnnotation):
    requires = DownloadDemoData
    envs = {
        "tool": "celltypist",
        "celltypist_args": {
            "model": MODEL, "over_clustering": False, "majority_voting": False
        },
    }


class CellTypeAnnotationSeurat(CellTypeAnnotation):
    requires = ModifyCellType
    envs = {
        "tool": "celltypist",
        "merge": True,
        "celltypist_args": {"model": MODEL},
    }


class Seurat2AnnData(Seurat2AnnData_):
    requires = CellTypeAnnotationSeurat


def pipeline():
    return get_pipeline(__file__).set_starts(DownloadDemoData)


def testing(pipen):
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "demo_2000_cells.annotated.h5ad",
        )
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
