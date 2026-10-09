"""Test SplitPipeSummary

This test is not run in CI (`LOCAL_ONLY`), since it requires `split-pipe` and
the sequencing data of two sublibraries. Put a small reference and two small
sublibraries in `tests/test_splitpipe/data/`:

- `ref.fa`, `ref.gtf`, `mkref.par`: the reference, used to build a genome with
    `SplitPipeRef` first
- `s1_R1.fastq.gz`, `s1_R2.fastq.gz`: the fastq files of the first sublibrary
- `s2_R1.fastq.gz`, `s2_R2.fastq.gz`: the fastq files of the second sublibrary

The data can be generated with
`python tests/test_splitpipe/data/make_test_data.py`. The barcodes in R2 are
drawn from the chemistry v2 whitelist, hence `chemistry = "v2"` below.

and run the test from the repo root:

    python tests/test_splitpipe/SplitPipeSummary/test.py
"""
from pathlib import Path

from biopipen.ns.splitpipe import (
    SplitPipeCombine,
    SplitPipeCount,
    SplitPipeRef,
    SplitPipeSummary,
)
from biopipen.core.testing import get_pipeline

DATA_DIR = Path(__file__).parent.parent / "data"


class SplitPipeRefTest(SplitPipeRef):
    input_data = [(str(DATA_DIR / "ref.fa"), str(DATA_DIR / "ref.gtf"))]
    # `parfile` keeps STAR's `genomeSAindexNbases` small for the tiny genome
    envs = {
        "genome_name": "test",
        "nthreads": 4,
        "parfile": str(DATA_DIR / "mkref.par"),
    }


class SplitPipeCountTest(SplitPipeCount):
    input = "fq1:files, fq2:files, id, genome_dir:dir"
    requires = SplitPipeRefTest
    input_data = lambda ch: [
        (
            [str(DATA_DIR / "s1_R1.fastq.gz")],
            [str(DATA_DIR / "s1_R2.fastq.gz")],
            "s1",
            ch.iloc[0, 0],
        ),
        (
            [str(DATA_DIR / "s2_R1.fastq.gz")],
            [str(DATA_DIR / "s2_R2.fastq.gz")],
            "s2",
            ch.iloc[0, 0],
        ),
    ]
    envs = {
        "chemistry": "v2",
        "kit": "WT",
        "nthreads": 4,
        "parfile": str(DATA_DIR / "count.par"),
    }
    forks = 2


class SplitPipeCombineTest(SplitPipeCombine):
    requires = SplitPipeCountTest
    input_data = lambda ch: [(list(ch.iloc[:, 0]), None)]
    envs = {
        "nthreads": 4,
    }


class SplitPipeSummaryTest(SplitPipeSummary):
    requires = SplitPipeCountTest


def pipeline():
    return get_pipeline(__file__, enable_report=True).set_start(SplitPipeRefTest)


def testing(pipen):
    outdir = Path(pipen.procs[-1].output_data.iloc[0, 0])
    metrics = outdir / "metrics.txt"
    assert metrics.is_file(), f"Metrics file was not created: {metrics}"
    assert metrics.stat().st_size > 0, f"Metrics file is empty: {metrics}"


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
