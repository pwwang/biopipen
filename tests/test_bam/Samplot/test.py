from pipen.channel import collapse_files, expand_dir
from biopipen.ns.web import Download
from biopipen.ns.bam import Samplot as Samplot_
from biopipen.ns.misc import Glob2Dir as Glob2Dir_
from biopipen.core.testing import get_pipeline

BAM_URLS = [
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12878_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12878_restricted.bam.bai",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12889_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12889_restricted.bam.bai",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12890_restricted.bam",
    "https://github.com/ryanlayer/samplot/raw/refs/heads/master/test/data/NA12890_restricted.bam.bai",
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
        "chrome": "chr4",
        "start": 115928726,
        "end": 115931880,
        # "sv_type": "DEL",
        "coverage_only": True,
    }


def pipeline():
    return get_pipeline(__file__).set_start(BamDownload)


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
