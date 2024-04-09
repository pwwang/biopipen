from pathlib import Path
from biopipen.core.testing import get_pipeline
from biopipen.ns.plot import ROC as _ROC


class ROC1(_ROC):
    envs = {"noids": True}


class ROC2(_ROC):
    envs = {"pos_label": "Ill"}


def pipeline():
    return (
        get_pipeline(__file__)
        .set_starts(ROC1, ROC2)
        .set_data(
            [Path(__file__).parent.joinpath("data", "single.txt")],
            [Path(__file__).parent.joinpath("data", "multi.txt")],
        )
    )


def testing(pipen):
    assert pipen._succeeded


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
