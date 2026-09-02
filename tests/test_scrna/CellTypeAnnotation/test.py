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
        "hitype": {"db": "hitypedb_pbmc3k"},
    }


class CellTypeAnnotationCell(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "cell",
        "set_ident": False,
        "cell_types": str(Path(__file__).parent / "data/celltype_annotation.tsv#1,2,3"),
    }


class CellTypeAnnotationScSorter(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "scsorter",
        "scsorter": {
            "db": str(
                Path(__file__).parent.parent
                / "Seurat"
                / "data/tcell.sccatch.RDS#celltype,gene"
            ),
        },
    }


class CellTypeAnnotationScType(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "sctype",
        "sctype": {
            "tissue": "Immune system",
            "db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
        },
    }


class CellTypeAnnotationSCINAUniversal(CellTypeAnnotation_):
    """SCINA with a universal marker table"""

    requires = PrepData
    envs = {
        "tool": "scina",
        "ident": "seurat_clusters",
        "scina": {"db": str(Path(__file__).parent / "data/markers.tsv")},
    }


class CellTypeAnnotationScSorterUniversal(CellTypeAnnotation_):
    """scSorter with a universal marker table + # column override"""

    requires = PrepData
    envs = {
        "tool": "scsorter",
        "scsorter": {
            "db": str(Path(__file__).parent / "data/markers.tsv#cell_type,gene"),
        },
    }


class CellTypeAnnotationScTypeUniversal(CellTypeAnnotation_):
    """sctype with a universal marker table (converted to ScType format)"""

    requires = PrepData
    envs = {
        "tool": "sctype",
        "sctype": {"db": str(Path(__file__).parent / "data/markers.tsv")},
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
        "ident": "seurat_clusters",
        "scina": {"db": str(Path(__file__).parent / "data/scina_signatures.csv")},
    }


class CellTypeAnnotationSingleR(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "singler",
        "singler": {"db": str(Path(__file__).parent / "data/singler_ref.rds")},
    }


class CellTypeAnnotationCelliD(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "cellid",
        "ident": "seurat_clusters",
        "cellid": {"db": str(Path(__file__).parent / "data/cellid_markers.csv")},
    }


class CellTypeAnnotationDeprecated(CellTypeAnnotation_):
    """Old-style flat envs must still work, with a deprecation warning"""

    requires = PrepData
    envs = {
        "tool": "hitype",
        "hitype_tissue": None,
        "hitype_db": "hitypedb_pbmc3k",
        "newcol": "CellType_old",
        "set_ident": False,
    }


class CellTypeAnnotationMultiCase(CellTypeAnnotation_):
    """Test multi-case annotation with two cluster-based tools"""

    requires = PrepData
    envs = {
        "cases": {
            "Sctype": {
                "tool": "sctype",
                "sctype": {
                    "tissue": "Immune system",
                    "db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
                },
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
                "sctype": {
                    "tissue": "Immune system",
                    "db": str(Path(__file__).parent / "data/ScTypeDB_short.xlsx"),
                },
                "anno_col": "CellType_sctype",
            },
            "Direct": {
                "tool": "direct",
                "cell_types": {
                    "c1": "Naive CD4+ T",
                    "c2": "NA",
                    "c5": "B",
                },
                "anno_col": "CellType_direct",
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
                "hitype": {"db": "hitypedb_pbmc3k"},
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


def testing(pipen):
    # assert pipen._succeeded

    # hitype (cluster-level): CellType column, no cell2celltype.tsv
    proc = get_proc(pipen, "CellTypeAnnotation")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotation")
    assert "CellType" in cols
    cluster_tsv = outfile.with_name(outfile.name + ".cluster2celltype.tsv")
    assert cluster_tsv.is_file()
    lines = cluster_tsv.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["Cluster", "Size", "DEFAULT"]
    assert sum(int(line.split("\t")[1]) for line in lines[1:]) == ncells
    assert not outfile.with_name(outfile.name + ".cell2celltype.tsv").exists()
    assert "SETTING IDENTS to 'CellType'" in proc.workdir.joinpath("0", "job.stdout").read_text()
    assert_idents_equal(pipen, "CellTypeAnnotation", "CellType")

    # cell tool with set_ident=False: original ident kept
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationCell")
    assert "CellType" in cols
    assert_idents_equal(pipen, "CellTypeAnnotationCell", "seurat_clusters")

    # SCINA (cell-level with ident): both cell and cluster outputs
    proc = get_proc(pipen, "CellTypeAnnotationSCINA")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationSCINA")
    assert "CellType" in cols
    assert "scina_celltype" in cols
    assert outfile.with_name(outfile.name + ".cluster2celltype.tsv").is_file()
    cell_tsv = outfile.with_name(outfile.name + ".cell2celltype.tsv")
    assert cell_tsv.is_file()
    lines = cell_tsv.read_text().splitlines()
    assert lines[0] == "Cell\tDEFAULT"
    assert len(lines) - 1 == ncells

    # CelliD (cell-level with ident): both cell and cluster outputs
    proc = get_proc(pipen, "CellTypeAnnotationCelliD")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationCelliD")
    assert "CellType" in cols
    assert "cellid_celltype" in cols
    assert outfile.with_name(outfile.name + ".cluster2celltype.tsv").is_file()
    cell_tsv = outfile.with_name(outfile.name + ".cell2celltype.tsv")
    assert cell_tsv.is_file()
    lines = cell_tsv.read_text().splitlines()
    assert lines[0] == "Cell\tDEFAULT"
    assert len(lines) - 1 == ncells
    # Values must be real cell types from the marker db, not barcodes.
    # Regression: RunCellHGT returns cell-types x cells, but annotate_cellid
    # used to assume cells x cell-types, yielding per-pathway garbage that
    # cbind recycled silently when ncells % ncelltypes == 0.
    cell_types = {"B cells", "DC", "Monocytes", "NK cells", "T cells"}
    assert set(line.split("\t")[1] for line in lines[1:]) <= cell_types
    clust_lines = outfile.with_name(
        outfile.name + ".cluster2celltype.tsv"
    ).read_text().splitlines()
    assert clust_lines[0].startswith("Cluster\t")
    assert set(line.split("\t")[-1] for line in clust_lines[1:]) <= cell_types

    # Universal marker table: SCINA (cell-level with ident)
    proc = get_proc(pipen, "CellTypeAnnotationSCINAUniversal")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationSCINAUniversal")
    assert "CellType" in cols
    assert "scina_celltype" in cols
    assert outfile.with_name(outfile.name + ".cluster2celltype.tsv").is_file()
    cell_tsv = outfile.with_name(outfile.name + ".cell2celltype.tsv")
    assert cell_tsv.is_file()
    lines = cell_tsv.read_text().splitlines()
    assert lines[0] == "Cell\tDEFAULT"
    assert len(lines) - 1 == ncells

    # Universal marker table: scSorter with # column override
    proc = get_proc(pipen, "CellTypeAnnotationScSorterUniversal")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationScSorterUniversal")
    assert "CellType" in cols
    assert outfile.with_name(outfile.name + ".cluster2celltype.tsv").is_file()

    # Universal marker table: sctype (converted to ScType format)
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationScTypeUniversal")
    assert "CellType" in cols

    # Old-style flat envs: deprecation warnings, newcol -> anno_col
    proc = get_proc(pipen, "CellTypeAnnotationDeprecated")
    stdout = proc.workdir.joinpath("0", "job.stdout").read_text()
    assert "`envs.hitype_db` is deprecated" in stdout
    assert "`envs.newcol` is deprecated" in stdout
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationDeprecated")
    assert "CellType_old" in cols
    assert_idents_equal(pipen, "CellTypeAnnotationDeprecated", "seurat_clusters")

    # Multi-case with add_prefix=False: per-case anno_col, last set_ident wins
    proc = get_proc(pipen, "CellTypeAnnotationMultiCase2")
    stdout = proc.workdir.joinpath("0", "job.stdout").read_text()
    assert "Multiple cases have set_ident enabled" in stdout
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationMultiCase2")
    assert "CellType_sctype" in cols
    assert "CellType_direct" in cols
    assert_idents_equal(pipen, "CellTypeAnnotationMultiCase2", "CellType_direct")


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
