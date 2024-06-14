from pathlib import Path

from biopipen.ns.regulation import (
    MotifAffinityTest as MotifAffinityTest_,
)
from biopipen.core.testing import get_pipeline

motiffile2 = Path(__file__).parent.parent / "MotifScan" / "data" / "motif2.txt"
motifdb = Path(__file__).parent.parent / "MotifScan" / "data" / "motifdb.txt"
reffa = Path(__file__).parent.parent.parent / "data" / "reference" / "hg19" / "chrs.fa"
varfile = Path(__file__).parent / "data" / "variants.bed"


class MotifAffinityTestMotifBreakR(MotifAffinityTest_):
    envs = {
        "tool": "motifbreakr",
        "regulator_col": 1,
        "motif_col": "Motif",
        "notfound": "ignore",
        "motifdb": motifdb,
    }


class MotifAffinityTestAtSNP(MotifAffinityTest_):
    envs = {
        "tool": "atsnp",
        "regulator_col": 1,
        "motif_col": "Motif",
        "notfound": "ignore",
        "motifdb": motifdb,
        "cutoff": 0.2,
        "atsnp_args": {"padj": "none"},
    }


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["+report"])
        get_pipeline(__file__)
        .set_starts(MotifAffinityTestMotifBreakR, MotifAffinityTestAtSNP)
        .set_data(
            [(motiffile2, varfile)],
            [(motiffile2, varfile)],
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
