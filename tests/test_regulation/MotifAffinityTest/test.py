from pathlib import Path

from biopipen.ns.regulation import (
    MotifAffinityTest as MotifAffinityTest_,
)
from biopipen.core.testing import get_pipeline

motiffile1 = Path(__file__).parent.parent / "MotifScan" / "data" / "motif1.txt"
motiffile2 = Path(__file__).parent.parent / "MotifScan" / "data" / "motif2.txt"
motifdb = Path(__file__).parent.parent / "MotifScan" / "data" / "motifdb.txt"
reffa = Path(__file__).parent.parent.parent / "data" / "reference" / "hg19" / "chrs.fa"
varfile = Path(__file__).parent / "data" / "variants.bed"


class MotifAffinityTestMotifBreakR(MotifAffinityTest_):
    envs = {
        "tool": "motifbreakr",
        "motif_col": "Motif",
        "notfound": "ignore",
        "motifdb": motifdb,
        "genome": "hg19",
    }


class MotifAffinityTestMotifBreakR_WithTF(MotifAffinityTest_):
    envs = {
        "tool": "motifbreakr",
        "motif_col": "Motif",
        "regulator_col": "TF",
        "notfound": "ignore",
        "motifdb": motifdb,
        "genome": "hg19",
    }


class MotifAffinityTestMotifBreakR_WithTFOnly(MotifAffinityTest_):
    envs = {
        "tool": "motifbreakr",
        "regulator_col": "TF",
        "notfound": "ignore",
        "motifdb": motifdb,
        "genome": "hg19",
        "regmotifs": motiffile1,
    }


class MotifAffinityTestAtSNP(MotifAffinityTest_):
    envs = {
        "tool": "atsnp",
        "motif_col": "Motif",
        "notfound": "ignore",
        "motifdb": motifdb,
        "cutoff": 0.2,
        "genome": "hg19",
        "atsnp_args": {"padj": "none"},
    }


class MotifAffinityTestAtSNP_WithTF(MotifAffinityTest_):
    envs = {
        "tool": "atsnp",
        "motif_col": "Motif",
        "regulator_col": "TF",
        "notfound": "ignore",
        "motifdb": motifdb,
        "cutoff": 0.2,
        "genome": "hg19",
        "atsnp_args": {"padj": "none"},
    }


class MotifAffinityTestAtSNP_WithTFOnly(MotifAffinityTest_):
    envs = {
        "tool": "atsnp",
        "regulator_col": "TF",
        "notfound": "ignore",
        "motifdb": motifdb,
        "cutoff": 0.2,
        "genome": "hg19",
        "atsnp_args": {"padj": "none"},
        "regmotifs": motiffile1,
    }


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["+report"])
        get_pipeline(__file__)
        .set_starts(
            MotifAffinityTestMotifBreakR,
            MotifAffinityTestMotifBreakR_WithTF,
            MotifAffinityTestMotifBreakR_WithTFOnly,
            MotifAffinityTestAtSNP,
            MotifAffinityTestAtSNP_WithTF,
            MotifAffinityTestAtSNP_WithTFOnly,
        )
        .set_data(
            [(motiffile2, varfile)],
            [(motiffile1, varfile)],
            [(motiffile1, varfile)],
            [(motiffile2, varfile)],
            [(motiffile1, varfile)],
            [(motiffile1, varfile)],
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
