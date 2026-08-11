from pipen.channel import collapse_files, expand_dir
from biopipen.ns.web import Download
from biopipen.ns.bam import (
    Samplot as Samplot_,
    BedtoolsCoverageBam as BedCoverageBam_,
    BedtoolsCoverageBamSummary as BedCoverageBamSummary_
)
from biopipen.ns.misc import Glob2Dir as Glob2Dir_, Str2File
from biopipen.core.testing import get_pipeline

BAM_URLS = [
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12878_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12878_restricted.bam.bai",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12889_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12889_restricted.bam.bai",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12890_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12890_restricted.bam.bai",
]


class BedFile(Str2File):
    input_data = [
        "chr4\t105928726\t115928725\n"
        "chr4\t115928726\t115931880\n"
    ]


class BamDownload(Download):
    input_data = BAM_URLS
    forks = 2


class Bams2Dir(Glob2Dir_):
    requires = BamDownload
    input_data = lambda ch: f"{collapse_files(ch).iloc[0, 0]}/*/output/*.bam*"


class Samplot(Samplot_):
    requires = Bams2Dir
    input_data = lambda ch: [list(expand_dir(ch, pattern="*.bam").iloc[:, 0])]
    envs = {
        "chrom": "chr4",
        "start": 115928726,
        "end": 115931880,
        # "sv_type": "DEL",
        "coverage_only": True,
    }


class BedCoverageBam(BedCoverageBam_):
    requires = Bams2Dir, BedFile
    input_data = lambda ch1, ch2: [
        (bam, ch2.iloc[0, 0])
        for bam in expand_dir(ch1, pattern="*.bam").iloc[:, 0]
    ]


class BedCoverageBamSummary(BedCoverageBamSummary_):
    requires = BedCoverageBam


def pipeline():
    return get_pipeline(__file__, forks=2).set_start(BedFile, BamDownload)


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
