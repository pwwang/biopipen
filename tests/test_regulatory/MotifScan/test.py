from pathlib import Path

from biopipen.ns.regulatory import (
    MotifScan as MotifScan_,
)
from biopipen.core.testing import get_pipeline

seqfile = Path(__file__).parent / "data" / "seq.fa"


class MotifScan1(MotifScan_):
    envs = {
        "motif_col": "Motif",
        "regulator_col": 1,
        "notfound": "ignore",
        "motifdb": Path(__file__).parent / "data" / "motifdb.txt",
        "cutoff": 0.64,
        "q": True,
        "q_cutoff": True,
    }


class MotifScan2(MotifScan_):
    envs = {
        "motifdb": Path(__file__).parent / "data" / "motifdb.txt",
        "cutoff": 5e-4,
    }


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["+report"])
        get_pipeline(__file__)
        .set_starts(MotifScan1, MotifScan2)
        .set_data(
            [(Path(__file__).parent / "data" / "motif1.txt", seqfile)],
            [(Path(__file__).parent / "data" / "motif2.txt", seqfile)],
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
