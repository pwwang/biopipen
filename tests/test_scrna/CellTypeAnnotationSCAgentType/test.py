"""Tests for CellTypeAnnotation with the `scagenttype` tool.

Local-only test (see run.env): it makes real LLM API calls via the
scAgentType agent, so it must not run in CI. It also requires the
`scagenttype` python package in the conda environment, installed from
source (not on PyPI), e.g.:
    pip install "scagenttype[llm] @ git+https://github.com/sathyasjali/scAgentType.git"

Run with an API key exported:
    OPENAI_API_KEY=sk-... python tests/test_scrna/CellTypeAnnotationSCAgentType/test.py
    # or via the test runner:
    OPENAI_API_KEY=sk-... bash tests/conda/run_test.sh \
        tests/test_scrna/CellTypeAnnotationSCAgentType FORCE=true

Copy the gitignored `.env` file from the CellTypeAnnotationLLM test dir (or
set the same values under the `OPENAI_*` names) to provide the API
credentials/base_url/model. This test reads `OPENAI_API_KEY`,
`OPENAI_BASE_URL`, and `OPENAI_MODEL` — the names scAgentType's OpenAI
backend reads natively. When copying from the LLM test dir, map
`LLMCELLTYPE_BASE_URL` -> `OPENAI_BASE_URL` and
`LLMCELLTYPE_MODEL` -> `OPENAI_MODEL`.

The dataset is a small pbmc3k subset (800 cells -> 8 clusters), so a single
agent run annotating 8 clusters is made. Cell-type names returned by the
model are not asserted — only the pipeline plumbing (CellType column,
cluster2celltype table, ident setting).
"""

import os
import subprocess

from pipen import Proc  # type: ignore
from biopipen.core.config import config
from biopipen.ns.scrna import CellTypeAnnotation as CellTypeAnnotation_
from biopipen.core.testing import get_pipeline
from dotenv import load_dotenv

env_path = os.path.join(os.path.dirname(__file__), ".env")
if os.path.isfile(env_path):
    load_dotenv(env_path)

BASE_URL = os.environ.get("OPENAI_BASE_URL")
API_KEY = os.environ.get("OPENAI_API_KEY")
if not API_KEY:
    raise SystemExit(
        "OPENAI_API_KEY is not set. This test makes real OpenAI API calls "
        "via the scAgentType python package and cannot run without it."
    )

MODEL = os.environ.get("OPENAI_MODEL", "gpt-4o-mini")


class PrepData(Proc):
    """Download pbmc3k and prepare a small clustered Seurat object"""

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
        # Keep it small: 800 cells -> 8 clusters (validated), so a single
        # scAgentType run on 8 clusters is made and the run stays quick.
        pbmc3k <- subset(pbmc3k, cells = sample(colnames(pbmc3k), 800))
        pbmc3k <- NormalizeData(pbmc3k)
        pbmc3k <- FindVariableFeatures(pbmc3k)
        pbmc3k <- ScaleData(pbmc3k)
        pbmc3k <- RunPCA(pbmc3k)
        pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
        pbmc3k <- FindClusters(pbmc3k, resolution = 0.8)
        saveRDS(pbmc3k, {{out.outfile | quote}})
    """


class CellTypeAnnotationSCAgentType(CellTypeAnnotation_):
    requires = PrepData
    envs = {
        "tool": "scagenttype",
        "scagenttype": {
            "api": "openai",
            "api_key": API_KEY,
            "model": MODEL,
            "base_url": BASE_URL,
            "python": "/home/pwwang/miniconda3/bin/python",
            "tissue": "Human peripheral blood",
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
        capture_output=True,
        text=True,
        check=True,
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
    # scagenttype (cluster-level): CellType column, no cell2celltype.tsv.
    # Cell-type names come from the model, so only assert the plumbing.
    proc = get_proc(pipen, "CellTypeAnnotationSCAgentType")
    outfile = proc.workdir.joinpath("0", "output", "pbmc3k.annotated.RDS")
    cols, idents, ncells = get_rds_info(pipen, "CellTypeAnnotationSCAgentType")
    assert "CellType" in cols
    cluster_tsv = outfile.with_name(
        outfile.name.replace(".RDS", "") + ".cluster2celltype.tsv"
    )
    assert cluster_tsv.is_file(), f"{cluster_tsv} does not exist"
    lines = cluster_tsv.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["Cluster", "Size", "DEFAULT"]
    assert sum(int(line.split("\t")[1]) for line in lines[1:]) == ncells
    assert all(line.split("\t")[2].strip() for line in lines[1:])
    assert not outfile.with_name(outfile.name + ".cell2celltype.tsv").exists()
    assert (
        "SETTING IDENTS to 'CellType'"
        in proc.workdir.joinpath("0", "job.stdout").read_text()
    )
    assert_idents_equal(pipen, "CellTypeAnnotationSCAgentType", "CellType")


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
