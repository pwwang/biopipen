from pathlib import Path

from biopipen.ns.misc import File2Proc
from biopipen.core.testing import get_pipeline


def pipeline():
    return get_pipeline(__file__).set_start(File2Proc).set_data([str(__file__)])

def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", Path(__file__).name)
    )
    assert outfile.exists()
    assert outfile.is_symlink()

if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
