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
        pbmc3k$Sample <- paste0("S", sample(1:12, nrow(pbmc3k), replace = TRUE))
        pbmc3k <- NormalizeData(pbmc3k)
        pbmc3k <- FindVariableFeatures(pbmc3k)
        pbmc3k <- ScaleData(pbmc3k)
        pbmc3k <- RunPCA(pbmc3k)
        pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
        pbmc3k <- FindClusters(pbmc3k, resolution = 0.8)
        pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)
        saveRDS(pbmc3k, {{out.outfile | quote}})
    """


class CellTypeAnnotation(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
    }


class CellTypeAnnotationCell(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "cell",
        "cell_types": str(Path(__file__).parent / "data/celltype_annotation.tsv#1,2,3"),
    }


class CellTypeAnnotationScSorter(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "scsorter",
        "scsorter_db": str(
            Path(__file__).parent.parent
            / "Seurat"
            / "data/tcell.sccatch.RDS#celltype,gene"
        ),
    }


class CellTypeAnnotationScType(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "sctype",
        "sctype_tissue": "Immune system",
        "sctype_db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
    }


class CellTypeAnnotationDirect(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "direct",
        "cell_types": [
            "Naive CD4+ T",
            "NA",
            "Memory CD4+",
            "Naive CD4+ T",
            "DC",
            "DC",
            "CD8+ T",
            "NK",
            "FCFR3A+ Mono",
            "CD8+ T",
        ],
        "merge": True,
    }


class CellTypeAnnotationDirect2(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "direct",
        "cell_types": {"c1": "Naive CD4+ T", "c2": "NA", "c5": "B"},
        "merge": True,
    }


class CellTypeAnnotationSCINA(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "scina",
        "scina_db": str(Path(__file__).parent / "data/scina_signatures.csv"),
    }


class CellTypeAnnotationSingleR(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "singler",
        "singler_db": str(Path(__file__).parent / "data/singler_ref.rds"),
    }


class CellTypeAnnotationCelliD(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "cellid",
        "cellid_db": str(Path(__file__).parent / "data/cellid_markers.csv"),
    }


class CellTypeAnnotationMultiCase(CellTypeAnnotation_):
    """Test multi-case annotation with two cluster-based tools"""

    requires = PrepData
    envs = {
        "cases": {
            "Sctype": {
                "tool": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
            },
            "Direct": {
                "tool": "direct",
                "cell_types": [
                    "Naive CD4+ T",
                    "NA",
                    "Memory CD4+",
                    "Naive CD4+ T",
                    "DC",
                    "DC",
                    "CD8+ T",
                    "NK",
                    "FCFR3A+ Mono",
                    "CD8+ T",
                ],
            },
        },
    }


class CellTypeAnnotationMultiCase2(CellTypeAnnotation_):
    """Test multi-case annotation with add_prefix=False"""

    requires = PrepData
    envs = {
        "add_prefix": False,
        "cases": {
            "Sctype": {
                "tool": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
                "newcol": "CellType_sctype",
            },
            "Direct": {
                "tool": "direct",
                "cell_types": {
                    "c1": "Naive CD4+ T",
                    "c2": "NA",
                    "c5": "B",
                },
                "newcol": "CellType_direct",
            },
        },
    }


class CellTypeAnnotationMultiCase3(CellTypeAnnotation_):
    """Test multi-case annotation with ncores"""

    requires = PrepData
    envs = {
        "ncores": 2,
        "cases": {
            "Hitype": {
                "tool": "hitype",
                "hitype_db": "hitypedb_pbmc3k",
            },
            "Direct": {
                "tool": "direct",
                "cell_types": [
                    "Naive CD4+ T",
                    "NA",
                    "Memory CD4+",
                    "Naive CD4+ T",
                    "DC",
                    "DC",
                    "CD8+ T",
                    "NK",
                    "FCFR3A+ Mono",
                    "CD8+ T",
                ],
            },
        },
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepData)


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
