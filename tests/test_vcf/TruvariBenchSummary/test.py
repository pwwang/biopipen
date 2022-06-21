from pathlib import Path

from datar.all import tibble
from pipen.channel import expand_dir
from biopipen.ns.vcf import TruvariBench, TruvariBenchSummary
from biopipen.ns.web import DownloadList
from biopipen.core.testing import get_pipeline


class DownloadList(DownloadList):
    """Download the vcfs"""


class TruvariBench(TruvariBench):
    requires = DownloadList
    input_data = lambda ch: tibble(
        compvcf=expand_dir(ch, pattern="input*.vcf.gz"),
        basevcf=expand_dir(ch, pattern="multi*.vcf.gz").iloc[0, 0],
    )
    envs = {
        "ref": str(
            Path(__file__).parent.parent.parent
            / "data"
            / "reference"
            / "hg19"
            / "chrs.fa"
        ),
    }


class TruvariBenchSummary(TruvariBenchSummary):
    requires = TruvariBench


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(DownloadList)
        .set_data([Path(__file__).parent / "data" / "vcfs.txt"])
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "truvari_bench.summary")
    )
    assert outfile.is_dir()

if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
