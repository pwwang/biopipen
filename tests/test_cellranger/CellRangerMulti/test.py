"""Test for CellRangerMulti and CellRangerMultiSummary.

Uses the same tiny nf-core test FASTQs as the CellRanger/test.py.
A "3' GEX only" configuration is described via envs (gex + libraries)
-- no manually-authored CSV needed. The multi config CSV is generated
automatically by the CellRangerMulti process.
The multi config references the same transcriptome as CellRangerCount via
config.ref.ref_cellranger_gex.
"""

from pathlib import Path
from biopipen.ns.cellranger import CellRangerMulti, CellRangerMultiSummary
from biopipen.ns.web import Download
from biopipen.core.config import config
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


class CellRangerMultiTest(CellRangerMulti):
    requires = DownloadData
    # Collect all 6 downloaded FASTQs into a single GEM well run.
    input_data = lambda ch: [(list(ch.iloc[:, 0]), "test_multi")]
    envs = {
        "ncores": 10,
        "cellranger": (
            "docker run -it --rm"
            f" -v {str(Path.cwd())}:{str(Path.cwd())}"
            " biopipen/cellranger:10.0.0"
        ),
        "gex": {
            "reference": config.ref.ref_cellranger_gex,
            "create_bam": False,
        },
        # Two GEX libraries from the same GEM well (Sample_X and Sample_Y)
        "libraries": [
            {"fastq_id": "Sample_X", "feature_types": "Gene Expression"},
            {"fastq_id": "Sample_Y", "feature_types": "Gene Expression"},
        ],
    }


class CellRangerMultiSummaryTest(CellRangerMultiSummary):
    requires = CellRangerMultiTest


def pipeline():
    return get_pipeline(__file__).set_start(DownloadData)


def testing(pipen): ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
