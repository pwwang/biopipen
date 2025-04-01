from biopipen.ns.web import Download
from biopipen.ns.cnv import (
    AneuploidyScore as AneuploidyScore_,
    AneuploidyScoreSummary as AneuploidyScoreSummary_,
    TMADScore as TMADScore_,
    TMADScoreSummary as TMADScoreSummary_,
)
from biopipen.core.testing import get_pipeline


CNS_URL = [
    "https://raw.githubusercontent.com/etal/cnvkit/refs/heads/master/test/formats/f-on-f.cns",  # noqa: E501
    "https://raw.githubusercontent.com/etal/cnvkit/refs/heads/master/test/formats/f-on-m.cns",  # noqa: E501
    "https://raw.githubusercontent.com/etal/cnvkit/refs/heads/master/test/formats/m-on-m.cns",  # noqa: E501
    "https://raw.githubusercontent.com/etal/cnvkit/refs/heads/master/test/formats/m-on-f.cns",  # noqa: E501
]

class AneuploidyScore(AneuploidyScore_):
    requires = Download
    forks = 4
    envs = {
        "cn_transform": "function(x) 2*2^x",
        "chrom_col": "chromosome",
        "start_col": "start",
        "end_col": "end",
        "seg_col": "log2",
        "genome": "hg19",
    }


class AneuploidyScoreSummary(AneuploidyScoreSummary_):
    requires = AneuploidyScore
    input_data = lambda ch: [(ch.outdir.tolist(), None)]


class TMADScore(TMADScore_):
    requires = Download
    forks = 4
    envs = {
        "segmean_transform": "function(x) 2*2^x",
        "chrom_col": "chromosome",
        "seg_col": "log2",
    }


class TMADScoreSummary(TMADScoreSummary_):
    requires = TMADScore
    input_data = lambda ch: [(ch.outfile.tolist(), None)]


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_start(Download)
        .set_data(CNS_URL)
    )


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
