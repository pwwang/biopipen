from pathlib import Path

from biopipen.ns.bed import BedConsensus
from biopipen.core.testing import get_pipeline

class BedConsensus(BedConsensus):

    envs = {
        "binsize": 400,
        "cutoff": 0.6,
        "genome": "hg19",
        "chrsize": Path(__file__).parent.parent.parent.joinpath(
            "data/reference/hg19/chrom.sizes"
        ).as_posix(),
    }

class BedConsensus1(BedConsensus):

    envs = {
        "binsize": 400,
        "cutoff": 0.9,
        "genome": "hg19",
        "ignore_scores": [0, 1, 2],
        "chrsize": Path(__file__).parent.parent.parent.joinpath(
            "data/reference/hg19/chrom.sizes"
        ).as_posix(),
    }


def pipeline():
    return get_pipeline(__file__, plugins=["no:report"]).set_starts(
        BedConsensus,
        BedConsensus1,
    ).set_data(
        [
            [
                Path(__file__).parent / "data" / "in1.bed",
                Path(__file__).parent / "data" / "in2.bed",
                Path(__file__).parent / "data" / "in3.bed",
            ]
        ]
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "in1_consensus.bed")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
