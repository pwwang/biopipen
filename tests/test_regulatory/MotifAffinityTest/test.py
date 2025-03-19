from pathlib import Path

from biopipen.ns.regulatory import (
    MotifAffinityTest as MotifAffinityTest_,
    VariantMotifPlot as VariantMotifPlot_,
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


class VariantMotifPlotBreakR(VariantMotifPlot_):
    requires = MotifAffinityTestMotifBreakR
    input_data = lambda ch: [ch.iloc[0, 0] / "motifbreakr.txt"]
    envs = {
        "genome": "hg19",
        "motifdb": motifdb,
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


class VariantMotifPlotBreakR_WithTF(VariantMotifPlot_):
    requires = MotifAffinityTestMotifBreakR_WithTF
    input_data = lambda ch: [ch.iloc[0, 0] / "motifbreakr.txt"]
    envs = {
        "genome": "hg19",
        "regulator_col": "Regulator",
        "motifdb": motifdb,
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


class VariantMotifPlotBreakR_WithTFOnly(VariantMotifPlot_):
    requires = MotifAffinityTestMotifBreakR_WithTFOnly
    input_data = lambda ch: [ch.iloc[0, 0] / "motifbreakr.txt"]
    envs = {
        "genome": "hg19",
        "motif_col": None,
        "regulator_col": "Regulator",
        "motifdb": motifdb,
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


class VariantMotifPlotAtSNP(VariantMotifPlot_):
    requires = MotifAffinityTestAtSNP
    input_data = lambda ch: [ch.iloc[0, 0] / "atsnp.txt"]
    envs = {
        "genome": "hg19",
        "motifdb": motifdb,
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


class VariantMotifPlotAtSNP_WithTF(VariantMotifPlot_):
    requires = MotifAffinityTestAtSNP_WithTF
    input_data = lambda ch: [ch.iloc[0, 0] / "atsnp.txt"]
    envs = {
        "genome": "hg19",
        "regulator_col": "Regulator",
        "motifdb": motifdb,
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


class VariantMotifPlotAtSNP_WithTFOnly(VariantMotifPlot_):
    requires = MotifAffinityTestAtSNP_WithTFOnly
    input_data = lambda ch: [ch.iloc[0, 0] / "atsnp.txt"]
    envs = {
        "genome": "hg19",
        "motif_col": None,
        "regulator_col": "Regulator",
        "motifdb": motifdb,
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
