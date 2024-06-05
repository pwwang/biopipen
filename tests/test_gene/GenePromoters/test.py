from pathlib import Path
from biopipen.ns.misc import Str2File
from biopipen.ns.gene import GenePromoters as GenePromoters_
from biopipen.core.testing import get_pipeline


class InputFile(Str2File):
    input_data = [
        (
            "Meta\tGene\n"
            "a\tCADPS2\n"
            "b\tABCA5\n"
            "c\tZNF451\n"
            "c1\tZNF451\n"
            "d\tPSAT1\n"
        )
    ]


class GenePromoters(GenePromoters_):
    requires = InputFile
    envs = {
        "genecol": 2,
        "refgene": Path(__file__).parent.parent.parent.joinpath(
            "data", "reference", "hg19", "refgene.gtf"
        )
    }


def pipeline():
    return get_pipeline(__file__).set_starts(InputFile)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
