from biopipen.ns.web import Download
from biopipen.ns.bam import CNVpytor
from biopipen.core.testing import get_pipeline


BAM_URL = (
    "https://github.com/VCCRI/SVPV/raw/master/example/NA12877_S1.partial.bam"
)


class CNVpytor(CNVpytor):
    requires = Download
    envs = {
        "binsizes": [10000],
    }


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(Download)
        .set_data([BAM_URL])
    )


def testing(pipen):
    assert pipen._succeeded


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
