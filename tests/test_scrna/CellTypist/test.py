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


class ModifyCellType(Proc):
    requires = DownloadDemoData
    input = "infile:file"
    output = "outfile:file:{{in.infile | stem}}.h5ad"
    lang = config.lang.python
    script = """
        import scanpy as sc
        infile = {{in.infile | r}}
        outfile = {{out.outfile | r}}
        obj = sc.read_h5ad(infile)
        obj.obs['cell_type'] = obj.obs['cell_type'].astype('category')
        old_cats = list(obj.obs['cell_type'].cat.categories)
        new_cats = [f"c{i+1}" for i in range(len(old_cats))]
        obj.obs['cell_type'] = obj.obs['cell_type'].cat.rename_categories(new_cats)
        obj.write_h5ad(outfile)
    """


class AnnData2Seurat(AnnData2Seurat_):
    requires = ModifyCellType
    envs = {"ident": "cell_type"}


class CellTypeAnnotationAnnData(CellTypeAnnotation):
    # No over_clustering, since no active_ident in the input
    # But majority_voting is default True
    # prediction will be in 'majority_voting' column
    requires = ModifyCellType
    envs = {
        "tool": "celltypist",
        "celltypist_args": {"model": MODEL},
    }


class CellTypeAnnotationAnnDataOverClustering(CellTypeAnnotation):
    requires = ModifyCellType
    envs = {
        "tool": "celltypist",
        "ident": "cell_type",
        "celltypist_args": {"model": MODEL},
    }


class CellTypeAnnotationSeurat(CellTypeAnnotation):
    requires = AnnData2Seurat
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

    # CellTypeAnnotationAnnData, because the input as no active_ident,
    # It is running no over_clustering by default
    cta_anndata_stdout = (
        [proc for proc in pipen.procs if proc.name == "CellTypeAnnotationAnnData"][0]
        .workdir.joinpath("0", "job.stdout")
    )
    assert "-c" not in cta_anndata_stdout.read_text()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
