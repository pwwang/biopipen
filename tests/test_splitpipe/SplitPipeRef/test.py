"""Test SplitPipeRef

This test is not run in CI (`LOCAL_ONLY`), since it requires `split-pipe` and
a reference genome. Put a small reference in `tests/test_splitpipe/data/`
(generated with `python tests/test_splitpipe/data/make_test_data.py`):

- `ref.fa`: the reference sequence in fasta format
- `ref.gtf`: the gene annotation in GTF format
- `mkref.par`: the parameter file with a small `genomeSAindexNbases` for the
    tiny genome

and run the test from the repo root:

    python tests/test_splitpipe/SplitPipeRef/test.py
"""
from pathlib import Path

from biopipen.ns.splitpipe import SplitPipeRef
from biopipen.core.testing import get_pipeline

DATA_DIR = Path(__file__).parent.parent / "data"


class SplitPipeRefTest(SplitPipeRef):
    input_data = [(str(DATA_DIR / "ref.fa"), str(DATA_DIR / "ref.gtf"))]
    # `mkref.par` sets `mkref_star_genomeSAindexNbases` for the tiny test
    # genome; STAR rejects its default of 14 for genomes this small.
    envs = {
        "genome_name": "test",
        "nthreads": 4,
        "parfile": str(DATA_DIR / "mkref.par"),
    }


def pipeline():
    return get_pipeline(__file__, enable_report=True).set_start(SplitPipeRefTest)


def testing(pipen):
    outdir = Path(pipen.procs[-1].output_data.iloc[0, 0])
    assert outdir.is_dir(), f"Genome directory was not created: {outdir}"
    assert any(outdir.iterdir()), f"Genome directory is empty: {outdir}"


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
