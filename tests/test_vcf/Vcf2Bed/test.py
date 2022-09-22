from pathlib import Path

from biopipen.ns.vcf import Vcf2Bed
from biopipen.core.testing import get_pipeline


INFILE = Path(__file__).parent.parent / "VcfFix" / "data" / "tofix.vcf"


def pipeline():
    return (
        get_pipeline(__file__, plugins=["no:report"])
        .set_start(Vcf2Bed)
        .set_data([str(INFILE)])
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "tofix.bed")
    )
    assert outfile.is_file()

if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
