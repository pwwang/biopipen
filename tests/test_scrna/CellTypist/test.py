from pathlib import Path

from biopipen.ns.scrna import CellTypeAnnotation, AnnData2Seurat, Seurat2AnnData
from biopipen.ns.web import Download
from biopipen.core.testing import get_pipeline

MODEL = Path(__file__).parent / "data" / "Immune_All_High.pkl"


class DownloadDemoData(Download):
    input_data = [
        "https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_2000_cells.h5ad"
    ]
    envs = {"tool": "aria2c"}


class AnnData2Seurat(AnnData2Seurat):
    requires = DownloadDemoData


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
    requires = AnnData2Seurat
    envs = {
        "tool": "celltypist",
        "celltypist_args": {"model": MODEL},
    }


class Seurat2AnnData(Seurat2AnnData):
    requires = CellTypeAnnotationSeurat


def pipeline():
    return get_pipeline(__file__).set_starts(DownloadDemoData)


def testing(pipen):
    assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "demo.2000.cells.annotated.h5ad",
        )
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
