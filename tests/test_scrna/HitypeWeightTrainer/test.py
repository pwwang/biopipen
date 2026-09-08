"""Tests for `HitypeWeightTrainer` — training hitype marker weights from a
Seurat object — and for the `CellTypeAnnotation` hitype tool consuming the
trained (weighted) universal marker tables.

Local-only test (see run.env): requires the `biopipen` R env with hitype >=
0.0.6 installed (`r-hitype`). Run locally with:

    python tests/test_scrna/HitypeWeightTrainer/test.py

Pipeline layout:

- `PrepData` prepares the pbmc3k data (SeuratData) with clustered Idents.
- `SyntheticData` creates a small two-type dataset (T_cell/B_cell, 100 cells
  each) whose cells express distinct marker gene sets (the same genes as
  `data/markers.tsv`) plus 40 noise genes, so the given-marker cases can be
  asserted semantically (labels T_cell/B_cell) instead of by cluster index.
- `HitypeWeightTrainerUniversal` trains from `data/markers.tsv` (a universal
  marker table) and `HitypeWeightTrainerNative` from `data/markers_native.txt`
  (the same markers in the hitype/ScType db format) — both must produce the
  same 9 marker rows with trained (not all 1) numeric weights.
- `HitypeWeightTrainerAuto` gives no `envs.markers`, so markers are discovered
  from pbmc3k by `hitype::find_markers()` over the `seurat_clusters` Idents,
  and weights are trained on them.
- `CellTypeAnnotationHitypeTrained` classifies the synthetic cells with the
  table trained in-pipeline by `HitypeWeightTrainerUniversal` (read from its
  deterministic job output via `envs.hitype.db`) — the **weighted** universal
  table through the fixed `CellTypeAnnotation` hitype route. This guards the
  regression where the `weight` column was dropped by the sctype-format
  conversion (`markers_to_sctype_df()`).
- `CellTypeAnnotationHitypeUniversal` does the same with the weightless
  `data/markers.tsv`, locking the same route for plain universal tables.
"""

import subprocess
from pathlib import Path

from pipen import Proc  # type: ignore
from biopipen.core.config import config
from biopipen.core.testing import _get_test_dirs, get_pipeline
from biopipen.ns.scrna import CellTypeAnnotation as CellTypeAnnotation_
from biopipen.ns.scrna import HitypeWeightTrainer as HitypeWeightTrainer_

HERE = Path(__file__).parent
MARKERS_TSV = HERE / "data/markers.tsv"
MARKERS_NATIVE_TXT = HERE / "data/markers_native.txt"
MARKERS_R = HERE.parents[2] / "biopipen" / "scripts" / "scrna" / "CellTypeAnnotation-markers.R"

# The given-markers trainer (HitypeWeightTrainerUniversal) writes its output
# at `<workdir>/<name>/HitypeWeightTrainerUniversal/0/output/syn.hitype.tsv`,
# which is where CellTypeAnnotationHitypeTrained reads it from
# (`envs.hitype.db`) after that trainer (a `requires` of the case) has run.
_NAME, _WORKDIR, _ = _get_test_dirs(__file__, False)
TRAINED_MARKERS = str(
    Path(_WORKDIR) / _NAME / "HitypeWeightTrainerUniversal" / "0" / "output"
    / "syn.hitype.tsv"
)

# Expected marker sets of the two synthetic cell types (data/markers.tsv and
# markers_native.txt are equivalent)
EXPECTED_MARKERS = {
    "T_cell": {"CD3D", "CD3E", "CD3G", "CD2", "IL7R"},
    "B_cell": {"MS4A1", "CD79A", "CD79B", "BANK1"},
}

CTA_LABELS = {"T_cell", "B_cell", "Unknown"}


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


class SyntheticData(Proc):
    """Two clean synthetic cell types (T_cell/B_cell, 100 cells each) whose
    cells express distinct marker gene sets at a high rate over a low
    background of 40 noise genes"""

    input = "name"
    input_data = ["syn"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        set.seed(8525)
        suppressMessages(library(Seurat))
        tmarkers <- c("CD3D", "CD3E", "CD3G", "CD2", "IL7R")
        bmarkers <- c("MS4A1", "CD79A", "CD79B", "BANK1")
        genes <- c(tmarkers, bmarkers, paste0("gene", 1:40))
        mkmat <- function(markers, n) {
            counts <- matrix(
                rpois(length(genes) * n, 0.3),
                length(genes), n,
                dimnames = list(genes, NULL)
            )
            counts[markers, ] <- counts[markers, ] +
                matrix(rpois(length(markers) * n, 12), length(markers))
            counts
        }
        counts <- cbind(mkmat(tmarkers, 100), mkmat(bmarkers, 100))
        colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
        sobj <- CreateSeuratObject(counts, min.cells = 1, min.features = 1)
        sobj$type <- rep(c("T_cell", "B_cell"), each = 100)
        sobj <- NormalizeData(sobj)
        # Cluster-level tools (`CellTypeAnnotation`'s hitype) group cells by
        # the ACTIVE Idents, and key their per-cluster mapping by the
        # `envs.ident` column — so both must carry the same values (as they
        # do for pbmc3k where Idents = seurat_clusters = `envs.ident`).
        Idents(sobj) <- sobj$type
        saveRDS(sobj, {{out.outfile | quote}})
    """


class HitypeWeightTrainerUniversal(HitypeWeightTrainer_):
    """Train from the universal marker table"""

    requires = SyntheticData
    envs = {
        "markers": str(MARKERS_TSV),
        "ident": "type",
    }


class HitypeWeightTrainerNative(HitypeWeightTrainer_):
    """Train from the same markers in the native hitype/ScType db format"""

    requires = SyntheticData
    envs = {
        "markers": str(MARKERS_NATIVE_TXT),
        "ident": "type",
    }


class HitypeWeightTrainerAuto(HitypeWeightTrainer_):
    """No `envs.markers`: find markers from the data with
    `hitype::find_markers()` over the `seurat_clusters` Idents"""

    requires = PrepData


class CellTypeAnnotationHitypeTrained(CellTypeAnnotation_):
    """Classify with the weighted table trained in-pipeline by
    `HitypeWeightTrainerUniversal`

    `requires` both SyntheticData (for the srtobj input) and
    HitypeWeightTrainerUniversal (so the weighted table below exists before
    this job runs; `srtobj` maps to SyntheticData's outfile positionally,
    the trainer's outfile column is unused).
    """

    requires = (SyntheticData, HitypeWeightTrainerUniversal)
    envs = {
        "tool": "hitype",
        "ident": "type",
        "hitype": {"db": TRAINED_MARKERS},
    }


class CellTypeAnnotationHitypeUniversal(CellTypeAnnotation_):
    """Classify with the weightless universal marker table (fixed canonical
    route for plain universal tables)"""

    requires = SyntheticData
    envs = {
        "tool": "hitype",
        "ident": "type",
        "hitype": {"db": str(MARKERS_TSV)},
    }


def pipeline():
    return get_pipeline(__file__).set_starts([PrepData, SyntheticData])


def get_proc(pipen, name):
    return [proc for proc in pipen.procs if proc.name == name][0]


def read_table(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return header, rows


def check_trainer(pipen, procname):
    """Common assertions for a given-markers trainer run: the output TSV has
    the universal-format header, exactly the fixture marker rows (one per
    gene), positive directions and trained numeric weights in [1, 5]."""
    proc = get_proc(pipen, procname)
    tsv = proc.workdir.joinpath("0", "output", "syn.hitype.tsv")
    assert tsv.is_file(), f"Missing trained marker table: {tsv}"
    header, rows = read_table(tsv)
    assert header == ["cell_type", "gene", "direction", "weight", "level"]
    pairs = {(row["cell_type"], row["gene"]) for row in rows}
    expected = {
        (ct, gene)
        for ct, genes in EXPECTED_MARKERS.items()
        for gene in genes
    }
    assert pairs == expected, f"Unexpected marker rows: {pairs ^ expected}"
    weights = [float(row["weight"]) for row in rows]
    # Some weights must actually have been learned (not all 1 — the
    # untrained fallback); calibrate against the real run if this bites
    assert len(set(round(w, 3) for w in weights)) >= 2, weights
    assert all(1 - 1e-6 <= w <= 5 + 1e-6 for w in weights), weights
    assert all(row["direction"] == "positive" for row in rows)
    assert all(row["level"] == "1" for row in rows)
    return tsv


def check_trainer_auto(pipen, procname):
    """The no-markers trainer: markers found by hitype::find_markers() over
    the pbmc3k clusters, then weighted — one row per found marker."""
    proc = get_proc(pipen, procname)
    tsv = proc.workdir.joinpath("0", "output", "pbmc3k.hitype.tsv")
    assert tsv.is_file(), f"Missing trained marker table: {tsv}"
    header, rows = read_table(tsv)
    assert header == ["cell_type", "gene", "direction", "weight", "level"]
    assert len(rows) > 0
    types = {row["cell_type"] for row in rows}
    # found markers come from the clusters (0..n) of the pbmc3k data
    assert len(types) >= 5, types
    weights = [float(row["weight"]) for row in rows]
    assert all(1 - 1e-6 <= w <= 5 + 1e-6 for w in weights), weights
    assert all(row["level"] == "1" for row in rows)
    return tsv


def get_rds_info(pipen, procname):
    """Return (meta.data columns, unique Idents, n cells) of the annotated RDS"""
    proc = get_proc(pipen, procname)
    rds = proc.workdir.joinpath("0", "output", "syn.annotated.RDS")
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


def check_cta_case(pipen, procname):
    """Common assertions for a CellTypeAnnotation (hitype, cluster-level) run
    on the synthetic data: annotated RDS with the CellType column (Idents set
    from it), a cluster2celltype.tsv of the right shape, and both synthetic
    cell types recovered (no all-Unknown / single-type collapse)."""
    proc = get_proc(pipen, procname)
    outprefix = proc.workdir.joinpath("0", "output", "syn.annotated")
    cols, idents, ncells = get_rds_info(pipen, procname)
    assert "CellType" in cols
    assert idents <= CTA_LABELS, f"Unexpected labels: {idents - CTA_LABELS}"
    # Labels of both synthetic types must be recovered (with clean synthetic
    # data, a wrong type assignment or an all-Unknown run is a regression)
    assert {"T_cell", "B_cell"} <= idents, idents
    cluster_tsv = outprefix.with_name(outprefix.name + ".cluster2celltype.tsv")
    assert cluster_tsv.is_file()
    assert not outprefix.with_name(outprefix.name + ".cell2celltype.tsv").exists()
    lines = cluster_tsv.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["Cluster", "Size", "DEFAULT"]
    assert sum(int(line.split("\t")[1]) for line in lines[1:]) == ncells
    script = f"""
        suppressMessages(library(Seurat))
        obj <- readRDS({str(proc.workdir.joinpath("0", "output", "syn.annotated.RDS"))!r})
        stopifnot(identical(
            as.character(Idents(obj)),
            as.character(obj@meta.data[["CellType"]])
        ))
    """
    subprocess.run(
        ["Rscript", "-e", script], capture_output=True, text=True, check=True
    )


def check_gs_weights_r():
    """The trained weights must survive into the hitype gene sets: run the
    same branch `CellTypeAnnotation`'s hitype route now takes for a universal
    marker table (version gate, apply_marker_filters, gs_prepare(df, NULL))
    on the trained TSV and assert the gene-set weights equal the file's
    `weight` column. A regression to markers_to_sctype_df() (which drops the
    column) would yield all-1 weights here."""
    script = f"""
        suppressMessages(library(biopipen.utils))
        suppressMessages(library(hitype))
        log <- get_logger()
        source({str(MARKERS_R)!r})
        stopifnot(packageVersion("hitype") >= "0.0.6")
        trained <- load_marker_table({str(TRAINED_MARKERS)!r})
        stopifnot(is_marker_canonical(trained))
        gs <- gs_prepare(trained, NULL)$gene_sets[[1]]
        stopifnot(!is.null(gs))
        # Every marker type must carry the file's numeric weights
        for (ct in names(gs)) {{
            rows <- trained[trained$cell_type == ct, , drop = FALSE]
            got <- gs[[ct]]$weights
            names(got) <- gs[[ct]]$markers
            exp <- setNames(as.numeric(rows$weight), rows$gene)
            stopifnot(identical(
                unname(got[names(exp)]),
                unname(exp)
            ))
        }}
        stopifnot(length(unique(unlist(lapply(gs, function(x) x$weights)))) >= 2)
        cat("GS_WEIGHTS_OK\\n")
    """
    result = subprocess.run(
        ["Rscript", "-e", script], capture_output=True, text=True, check=True
    )
    assert "GS_WEIGHTS_OK" in result.stdout


# pyright: reportOperatorIssue=false
def testing(pipen):
    check_trainer(pipen, "HitypeWeightTrainerUniversal")
    check_trainer(pipen, "HitypeWeightTrainerNative")
    check_trainer_auto(pipen, "HitypeWeightTrainerAuto")

    check_cta_case(pipen, "CellTypeAnnotationHitypeTrained")
    check_cta_case(pipen, "CellTypeAnnotationHitypeUniversal")

    check_gs_weights_r()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
