from biopipen.ns.scrna import (
    LoomTo10X as LoomTo10X_,
    SeuratPreparing as SeuratPreparing_,
)
from biopipen.core.proc import Proc
from biopipen.ns.web import Download
from biopipen.core.testing import get_pipeline


class DownloadLoom(Download):
    """Download the loom file"""

    input_data = [
        "http://loom.linnarssonlab.org/clone/Mousebrain.org.level6/"
        "L6_Peripheral_sensory_peptidergic_neurons.loom"
    ]


class LoomTo10X(LoomTo10X_):
    requires = DownloadLoom


class CreateMetaFile(Proc):
    """Create meta file for SeuratPreparing"""
    requires = LoomTo10X
    input = "indir:dir"
    output = "outfile:file:meta.txt"
    script = """
        indir={{in.indir | quote}}
        outfile={{out.outfile | quote}}
        echo -e "Sample\\tRNAData" > $outfile
        echo -e "Sample1\\t$indir" >> $outfile
    """


# check if the conversion is correct
class SeuratPreparing(SeuratPreparing_):
    requires = CreateMetaFile
    # It's just a single sample
    envs = {"no_integration": True}


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_starts(DownloadLoom)
    )


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
