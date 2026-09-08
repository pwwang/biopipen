"""Tests for CellTypeAnnotation with the `garnett` tool.

Local-only test (see run.env): `garnett::classify_cells()` needs the
`cole-trapnell-lab/garnett` package (monocle3 branch) plus `monocle3` and
`SeuratWrappers`, which are not in the CI environment. They are installed ad
hoc into the local conda env (`conda install -c bioconda r-monocle3`, then
`remotes::install_github("cole-trapnell-lab/garnett", ref = "monocle3")` and
`satijalab/SeuratWrappers`). Run locally with:

    python tests/test_scrna/CellTypeAnnotationGarnett/test.py
    # or via the test runner:
    bash tests/conda/run_test.sh tests/test_scrna/CellTypeAnnotationGarnett FORCE=true

The classifier is the official pre-trained human PBMC classifier from
<https://cole-trapnell-lab.github.io/garnett/classifiers/> (trained on ENSEMBL
gene IDs), so `db: org.Hs.eg.db` + `cds_gene_id_type: SYMBOL` convert the
gene-symbol expression data into the classifier's ID space. Cell-type labels
are asserted against the classifier's known label universe, with at least one
cell/cluster annotated (not all `Unknown`), which catches wrong gene-ID
conversion regressions.
"""

import subprocess
from pathlib import Path

from pipen import Proc  # type: ignore
from biopipen.core.config import config
from biopipen.ns.scrna import CellTypeAnnotation as CellTypeAnnotation_
from biopipen.core.testing import get_pipeline


class PrepData(Proc):
    """Download pbmc3k, normalize, and cluster"""

    input = "name"
    input_data = ["pbmc3k"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        set.seed(8525)
        options(timeout = 600)
        library(Seurat)
        library(SeuratData)
        tryCatch({
            InstallData("pbmc3k")
        }, error = function(e) {
            # https://github.com/satijalab/seurat-data/issues/23#issuecomment-1227111059
            install.packages(
                "pbmc3k.SeuratData",
                repos = "http://seurat.nygenome.org/",
                type = "source"
            )
        })
        pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
        pbmc3k <- NormalizeData(pbmc3k)
        pbmc3k <- FindVariableFeatures(pbmc3k)
        pbmc3k <- ScaleData(pbmc3k)
        pbmc3k <- RunPCA(pbmc3k)
        pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
        pbmc3k <- FindClusters(pbmc3k, resolution = 0.8)
        saveRDS(pbmc3k, {{out.outfile | quote}})
    """


class CellTypeAnnotationGarnett(CellTypeAnnotation_):
    """Garnett with the official hsPBMC classifier (ENSEMBL ids -> symbols)"""

    requires = PrepData
    envs = {
        "tool": "garnett",
        "ident": "seurat_clusters",
        "garnett": {
            "classifier": str(Path(__file__).parent / "data/hsPBMC_20191017.RDS"),
            "db": "org.Hs.eg.db",
            "cds_gene_id_type": "SYMBOL",
        },
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepData)


def get_proc(pipen, name):
    return [proc for proc in pipen.procs if proc.name == name][0]


def get_rds_info(pipen, procname):
    """Return (meta.data columns, unique Idents, n cells) of the annotated RDS"""
    proc = get_proc(pipen, procname)
    rds = proc.workdir.joinpath("0", "output", "pbmc3k.annotated.RDS")
    script = f"""
        suppressMessages(library(Seurat))
        obj <- readRDS({str(rds)!r})
        cat("COLS:", paste(colnames(obj@meta.data), collapse = ","), "\\n", sep = "")
        cat("IDENTS:", paste(unique(as.character(Idents(obj))), collapse = ","), "\\n", sep = "")
        cat("NCELLS:", nrow(obj@meta.data), "\\n", sep = "")
    """
    result = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True, text=True, check=True,
    )
    cols = None
    idents = None
    ncells = None
    for line in result.stdout.splitlines():
        if line.startswith("COLS:"):
            cols = set(line[5:].split(","))
        elif line.startswith("IDENTS:"):
            idents = set(line[7:].split(","))
        elif line.startswith("NCELLS:"):
            ncells = int(line[7:])
    return cols, idents, ncells


def assert_idents_equal(pipen, procname, colname):
    """Assert the final Idents equal the values of meta.data[[colname]]"""
    proc = get_proc(pipen, procname)
    rds = proc.workdir.joinpath("0", "output", "pbmc3k.annotated.RDS")
    script = f"""
        suppressMessages(library(Seurat))
        obj <- readRDS({str(rds)!r})
        stopifnot(identical(
            as.character(Idents(obj)),
            as.character(obj@meta.data[[{colname!r}]])
        ))
    """
    subprocess.run(
        ["Rscript", "-e", script], capture_output=True, text=True, check=True
    )


# pyright: reportOperatorIssue=false
def testing(pipen):
    # garnett (cell-level with ident): both cell and cluster outputs
    proc = get_proc(pipen, "CellTypeAnnotationGarnett")
    outprefix = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationGarnett")
    assert "CellType" in cols
    assert "garnett_celltype" in cols
    cell_tsv = outprefix.with_name(outprefix.name + ".cell2celltype.tsv")
    cluster_tsv = outprefix.with_name(outprefix.name + ".cluster2celltype.tsv")
    assert cell_tsv.is_file()
    assert cluster_tsv.is_file()
    lines = cell_tsv.read_text().splitlines()
    assert lines[0] == "Cell\tDEFAULT"
    assert len(lines) - 1 == ncells
    # Labels are from the hsPBMC classifier's known universe; some cells must
    # be annotated (not all Unknown) or the gene-ID conversion is broken
    labels = {
        "T cells", "CD4 T cells", "CD8 T cells", "B cells", "Monocytes",
        "NK cells", "Dendritic cells", "CD34+", "Unknown",
    }
    values = set(line.split("\t")[1] for line in lines[1:])
    assert values <= labels, f"Unexpected labels: {values - labels}"
    assert any(v != "Unknown" for v in values)
    clust_lines = cluster_tsv.read_text().splitlines()
    assert clust_lines[0].startswith("Cluster\t")
    cl_values = set(line.split("\t")[-1] for line in clust_lines[1:])
    assert cl_values <= labels, f"Unexpected cluster labels: {cl_values - labels}"
    assert any(v != "Unknown" for v in cl_values)
    assert (
        "SETTING IDENTS to 'CellType'"
        in proc.workdir.joinpath("0", "job.stdout").read_text()
    )
    assert_idents_equal(pipen, "CellTypeAnnotationGarnett", "CellType")


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
