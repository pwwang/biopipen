from pathlib import Path

from biopipen.ns.tcr import TESSA as TESSA_
from biopipen.core.testing import get_pipeline


class TESSA(TESSA_):
    envs = {
        "prefix": "",
        "max_iter": 100,
    }


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["no:report"])
        get_pipeline(__file__)
        .set_starts(TESSA)
        .set_data(
            [
                (
                    Path(__file__).parent / "data" / "tcrdata.RDS",
                    Path(__file__).parent / "data" / "expdata.txt.gz",
                )
            ]
        )
    )


def testing(pipen):
    assert pipen._succeeded
    outfile = pipen.procs[-1].workdir.joinpath(
        "0",
        "output",
        "tcrdata.tessa.txt",
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
