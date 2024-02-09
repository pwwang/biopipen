from pathlib import Path

from biopipen.ns.rnaseq import UnitConversion as UnitConversion_
from biopipen.core.testing import get_pipeline


class UnitConversion(UnitConversion_):
    envs = {
        "inunit": "count",
        "outunit": "cpm",
    }


class UnitConversion2(UnitConversion_):
    envs = {
        "inunit": "log2(count + 1)",
        "outunit": "log2(cpm + 1)",
    }


class UnitConversion3(UnitConversion_):
    envs = {
        "inunit": "count",
        "outunit": "log2(count + 1)",
    }


def pipeline():
    return (
        get_pipeline(__file__, plugins=["no:report"])
        # get_pipeline(__file__)
        .set_starts(UnitConversion, UnitConversion2, UnitConversion3)
        .set_data(
            [Path(__file__).parent / "data" / "exprs.txt"],
            [Path(__file__).parent / "data" / "exprs.txt"],
            [Path(__file__).parent / "data" / "exprs.txt"],
        )
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "exprs.txt")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
