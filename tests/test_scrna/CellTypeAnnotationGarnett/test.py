"""Tests for CellTypeAnnotation with the `garnett` tool, and for the
`GarnettClassifierTrainer` process that trains the classifiers it consumes.

Local-only test (see run.env): `garnett` (monocle3 branch), `monocle3` and
`SeuratWrappers` are not in the CI environment. They are installed ad hoc
into the local conda env (`conda install -c bioconda r-monocle3`, then
`remotes::install_github("cole-trapnell-lab/garnett", ref = "monocle3")` and
`satijalab/SeuratWrappers`). Run locally with:

    python tests/test_scrna/CellTypeAnnotationGarnett/test.py
    # or via the test runner:
    bash tests/conda/run_test.sh tests/test_scrna/CellTypeAnnotationGarnett FORCE=true

Pipeline layout (all from the same pbmc3k `PrepData`):

- `GarnettClassifierTrainer` trains from a universal marker table
  (auto-converted to a garnett marker file; `direction: negative` markers
  becoming `not expressed:` rules) and `GarnettClassifierTrainerNative` from
  an equivalent garnett-native marker file (used as is — no converted side
  file, proving the format dispatch by content sniffing).
- `CellTypeAnnotationGarnett` classifies with the official pre-trained human
  PBMC classifier from <https://cole-trapnell-lab.github.io/garnett/classifiers/>
  (trained on ENSEMBL gene IDs), so `db: org.Hs.eg.db` +
  `cds_gene_id_type: SYMBOL` convert the gene-symbol expression data into the
  classifier's ID space. Cell-type labels are asserted against the
  classifier's known label universe, with at least one cell/cluster annotated
  (not all `Unknown`), which catches wrong gene-ID conversion regressions.
- `CellTypeAnnotationGarnettTrained` classifies with the classifier trained
  by `GarnettClassifierTrainer` in this same pipeline (a `requires` of the
  case), reading it from its deterministic job output via `envs.garnett.classifier`.
  It is trained on the very same pbmc3k data, so clean markers make the
  cv.glmnet fit near-perfectly separable: lambda lands at the floor of
  garnett's grid, the huge coefficients overflow `exp()` in the softmax, and
  most cells come back `Unknown`. Only a non-empty classified set is
  asserted here — the all-`Unknown` glmnet >= 4.0 regression is still caught —
  while the official-classifier case above guards healthy classification.
"""

import subprocess
from pathlib import Path

from pipen import Proc  # type: ignore
from biopipen.core.config import config
from biopipen.ns.scrna import CellTypeAnnotation as CellTypeAnnotation_
from biopipen.ns.scrna import GarnettClassifierTrainer as GarnettClassifierTrainer_
from biopipen.core.testing import _get_test_dirs, get_pipeline

HERE = Path(__file__).parent
REPO = HERE.parents[2]
MARKERS_R = REPO / "biopipen" / "scripts" / "scrna" / "CellTypeAnnotation-markers.R"
UNIVERSAL_TSV = HERE / "data" / "garnett_markers.tsv"
NATIVE_TXT = HERE / "data" / "garnett_markers_native.txt"

# name/workdir of this test's pipeline, as `pipeline()` builds it below. The
# universal-trainer job writes its classifier at
# `<workdir>/<name>/GarnettClassifierTrainer/0/output/pbmc3k.classifier.RDS`,
# which is where `CellTypeAnnotationGarnettTrained` reads it from
# (`envs.garnett.classifier`), after that trainer (a `requires` of the case)
# has run.
_NAME, _WORKDIR, _ = _get_test_dirs(__file__, False)
TRAINED_CLASSIFIER = str(
    Path(_WORKDIR) / _NAME / "GarnettClassifierTrainer" / "0" / "output"
    / "pbmc3k.classifier.RDS"
)

# Expected content of the converted garnett marker file (side output of the
# universal-table input): one `> <cell type>` block per table cell type, in
# table order, genes sorted, negatives as `not expressed:`.
EXPECTED_MARKERS_TXT = "\n".join([
    "> T_cell",
    "expressed: CD3D, CD3E, CD3G",
    "> B_cell",
    "expressed: CD79A, MS4A1",
    "not expressed: CD3D",
    "> Monocyte",
    "expressed: CD14, LYZ",
    "> NK_cell",
    "expressed: KLRD1, NKG7",
])

# Label universe of the in-pipeline-trained classifiers (marker cell types
# plus garnett's `Unknown` sentinel)
TRAINED_LABELS = {"T_cell", "B_cell", "Monocyte", "NK_cell", "Unknown"}


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


class GarnettClassifierTrainer(GarnettClassifierTrainer_):
    """Train from the universal marker table (auto-converted)"""

    requires = PrepData
    input_data = lambda ch: [  # ch: PrepData's output data
        (ch.outfile[0], str(UNIVERSAL_TSV))
    ]
    envs = {
        "seed": 8525,
    }


class GarnettClassifierTrainerNative(GarnettClassifierTrainer_):
    """Train from the garnett-native marker file (used as is)"""

    requires = PrepData
    input_data = lambda ch: [  # ch: PrepData's output data
        (ch.outfile[0], str(NATIVE_TXT))
    ]
    envs = {
        "seed": 8525,
    }


class CellTypeAnnotationGarnett(CellTypeAnnotation_):
    """Garnett with the official hsPBMC classifier (ENSEMBL ids -> symbols)"""

    requires = PrepData
    envs = {
        "tool": "garnett",
        "ident": "seurat_clusters",
        "garnett": {
            "classifier": str(HERE / "data/hsPBMC_20191017.RDS"),
            "db": "org.Hs.eg.db",
            "cds_gene_id_type": "SYMBOL",
        },
    }


class CellTypeAnnotationGarnettTrained(CellTypeAnnotation_):
    """Garnett with the classifier trained in-pipeline by GarnettClassifierTrainer

    `requires` both PrepData (for the srtobj input) and GarnettClassifierTrainer
    (so the classifier below exists before this job runs; `srtobj` maps to
    PrepData's outfile positionally, the trainer's outfile column is unused).
    The classifier is trained with `db: "none"` (gene symbols kept as is), so
    classification needs no ID conversion either.
    """

    requires = (PrepData, GarnettClassifierTrainer)
    envs = {
        "tool": "garnett",
        "ident": "seurat_clusters",
        "garnett": {
            "classifier": TRAINED_CLASSIFIER,
            "db": "none",
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


def check_trainer(pipen, procname, expected_side_markers):
    """Check a trained classifier: RDS well-formed, tree has the 4 marker
    cell types, and the marker conversion side file (or its absence) matches
    the input format."""
    proc = get_proc(pipen, procname)
    outdir = proc.workdir.joinpath("0", "output")
    classifier = outdir.joinpath("pbmc3k.classifier.RDS")
    assert classifier.is_file(), f"Missing classifier: {classifier}"
    side_markers = outdir.joinpath("pbmc3k.classifier.markers.txt")
    assert side_markers.is_file() == expected_side_markers, (
        f"Unexpected converted marker side file: {side_markers}"
    )
    if expected_side_markers:
        assert side_markers.read_text() == EXPECTED_MARKERS_TXT + "\n"
    script = f"""
        suppressMessages(library(garnett))
        suppressMessages(library(biopipen.utils))
        log <- get_logger()
        classifier <- read_obj({str(classifier)!r})
        stopifnot(inherits(classifier, "garnett_classifier"))
        tree <- classifier@classification_tree
        vnames <- igraph::V(tree)$name
        cat("TREE:", paste(sort(vnames), collapse = ","), "\\n", sep = "")
        # at least one node was actually trained (has a fitted model)
        stopifnot(any(!vapply(
            igraph::vertex_attr(tree, "model"), is.null, logical(1)
        )))
    """
    result = subprocess.run(
        ["Rscript", "-e", script], capture_output=True, text=True, check=True
    )
    tree = None
    for line in result.stdout.splitlines():
        if line.startswith("TREE:"):
            tree = set(line[5:].split(","))
    assert tree >= {"T_cell", "B_cell", "Monocyte", "NK_cell"}, tree


def check_classifier_case(pipen, procname, labels, cell_annotated=True):
    """Common assertions for a CellTypeAnnotation (garnett) run: annotated RDS
    with both annotation columns, both cell/cluster TSVs of the right shape,
    labels within the classifier's universe, and Idents set from CellType."""
    proc = get_proc(pipen, procname)
    outprefix = proc.workdir.joinpath("0", "output", "pbmc3k.annotated")
    cols, idents, ncells = get_rds_info(pipen, procname)
    assert "CellType" in cols
    assert "garnett_celltype" in cols
    cell_tsv = outprefix.with_name(outprefix.name + ".cell2celltype.tsv")
    cluster_tsv = outprefix.with_name(outprefix.name + ".cluster2celltype.tsv")
    assert cell_tsv.is_file()
    assert cluster_tsv.is_file()
    lines = cell_tsv.read_text().splitlines()
    assert lines[0] == "Cell\tDEFAULT"
    assert len(lines) - 1 == ncells
    values = set(line.split("\t")[1] for line in lines[1:])
    assert values <= labels, f"Unexpected labels: {values - labels}"
    if cell_annotated:
        # Some cells must be annotated (not all Unknown) or the gene-ID
        # conversion / glmnet workaround is broken
        assert any(v != "Unknown" for v in values)
    clust_lines = cluster_tsv.read_text().splitlines()
    assert clust_lines[0].startswith("Cluster\t")
    cl_values = set(line.split("\t")[-1] for line in clust_lines[1:])
    assert cl_values <= labels, f"Unexpected cluster labels: {cl_values - labels}"
    assert (
        "SETTING IDENTS to 'CellType'"
        in proc.workdir.joinpath("0", "job.stdout").read_text()
    )
    assert_idents_equal(pipen, procname, "CellType")


def check_markers_r():
    """Unit-check the shared marker helpers of CellTypeAnnotation-markers.R."""
    script = f"""
        suppressMessages(library(biopipen.utils))
        log <- get_logger()
        source({str(MARKERS_R)!r})

        # is_garnett_native_marker(): content sniffing
        stopifnot(identical(
            is_garnett_native_marker({str(NATIVE_TXT)!r}), TRUE
        ))
        stopifnot(identical(
            is_garnett_native_marker({str(UNIVERSAL_TSV)!r}), FALSE
        ))
        stopifnot(identical(
            is_garnett_native_marker("/no/such/file.markers.txt"), FALSE
        ))
        rds <- tempfile(fileext = ".rds")
        saveRDS(data.frame(cell_type = "T_cell", gene = "CD3D"), rds)
        stopifnot(identical(is_garnett_native_marker(rds), FALSE))
        unlink(rds)

        # markers_to_garnett_file(): species filter + ordering
        tbl <- tempfile(fileext = ".tsv")
        writeLines(c(
            "cell_type\\tgene\\tspecies",
            "T_cell\\tCD3D\\tHuman",
            "T_cell\\tCD3E\\tMouse",
            "B_cell\\tMS4A1\\tMouse",
            "B_cell\\tMS4A1\\tHuman",
            "B_cell\\tCD79A\\tMouse"
        ), tbl)
        out <- tempfile(fileext = ".markers.txt")
        markers_to_garnett_file(load_marker_table(tbl), out, species = "Mouse")
        stopifnot(identical(readLines(out), c(
            "> T_cell",
            "expressed: CD3E",
            "> B_cell",
            "expressed: CD79A, MS4A1"
        )))
        unlink(c(tbl, out))

        # markers_to_garnett_file(): a cell type without positive markers errors
        tbl2 <- tempfile(fileext = ".tsv")
        writeLines(c(
            "cell_type\\tgene\\tdirection",
            "T_cell\\tCD3D\\tnegative"
        ), tbl2)
        out2 <- tempfile(fileext = ".markers.txt")
        err <- tryCatch(
            markers_to_garnett_file(load_marker_table(tbl2), out2),
            error = identity
        )
        stopifnot(
            inherits(err, "error"),
            grepl("no positive markers", conditionMessage(err))
        )
        unlink(c(tbl2, out2))
    """
    subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True, text=True, check=True,
    )


# pyright: reportOperatorIssue=false
def testing(pipen):
    check_markers_r()
    check_trainer(pipen, "GarnettClassifierTrainer", expected_side_markers=True)
    check_trainer(pipen, "GarnettClassifierTrainerNative", expected_side_markers=False)

    # Official classifier: healthy classification (labels from the hsPBMC
    # classifier's known universe, most cells annotated)
    check_classifier_case(
        pipen,
        "CellTypeAnnotationGarnett",
        {
            "T cells", "CD4 T cells", "CD8 T cells", "B cells", "Monocytes",
            "NK cells", "Dendritic cells", "CD34+", "Unknown",
        },
    )

    # In-pipeline-trained classifier: trained on the same data, so most cells
    # come back Unknown (see the module docstring) — guard the all-Unknown
    # glmnet >= 4.0 regression with at least one annotated cell
    check_classifier_case(
        pipen, "CellTypeAnnotationGarnettTrained", TRAINED_LABELS
    )


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
