from pathlib import Path

from biopipen.ns.vcf import Vcf2Bed
from biopipen.core.testing import get_pipeline


INFILE = Path(__file__).parent / "data" / "tobed.vcf"


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(Vcf2Bed)
        .set_data([str(INFILE)])
    )


def testing(pipen):
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "tobed.bed")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
