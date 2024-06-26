import sys
from pathlib import Path

from biopipen.ns.tcr import CloneSizeQQPlot
from biopipen.core.testing import get_pipeline

sys.path.insert(0, str(Path(__file__).parent.parent))
from TCRClusterStats.test import PrepareImmdata  # noqa: E402


class PrepareImmdata(PrepareImmdata):
    ...


class CloneSizeQQPlot(CloneSizeQQPlot):
    requires = PrepareImmdata
    envs = {
        "subject": "Status",
        "group": "ID",
        "order": ["C1", "C2", "MS1", "MS2"],
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareImmdata)


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
