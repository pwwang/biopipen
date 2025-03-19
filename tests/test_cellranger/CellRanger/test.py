from biopipen.ns.cellranger import CellRangerCount, CellRangerSummary
from biopipen.ns.web import Download
from biopipen.core.testing import get_pipeline


DATA_BASE_URL = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/"
    "scrnaseq/testdata/cellranger"
)
DATA_FILES = (
    "Sample_X_S1_L001_R1_001.fastq.gz",
    "Sample_X_S1_L001_R2_001.fastq.gz",
    "Sample_Y_S1_L001_R1_001.fastq.gz",
    "Sample_Y_S1_L001_R2_001.fastq.gz",
    "Sample_Y_S1_L002_R1_001.fastq.gz",
    "Sample_Y_S1_L002_R2_001.fastq.gz",
)


class DownloadData(Download):
    input_data = [f"{DATA_BASE_URL}/{file}" for file in DATA_FILES]
    envs = {"tool": "aria2c"}


class CellRangerCountTest(CellRangerCount):
    requires = DownloadData
    input_data = lambda ch: [
        (ch.iloc[:2, 0], "Sample_1"),
        (ch.iloc[2:, 0], "Sample_2"),
    ]
    envs = {"ncores": 10}
    forks = 2


class CellRangerSummaryTest(CellRangerSummary):
    requires = CellRangerCountTest


def pipeline():
    return get_pipeline(__file__).set_start(DownloadData)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
