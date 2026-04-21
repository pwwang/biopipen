"""Test for CellRangerMulti and CellRangerMultiSummary.

Uses the same tiny nf-core test FASTQs as the CellRanger/test.py.
A "3' GEX only" multi config CSV is generated dynamically from the
downloaded FASTQs (gene expression libraries only, no VDJ/Feature Barcode).
The multi config references the same transcriptome as CellRangerCount via
config.ref.ref_cellranger_gex.
"""

from biopipen.ns.cellranger import CellRangerMulti, CellRangerMultiSummary
from biopipen.ns.web import Download
from biopipen.core.proc import Proc
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


class CreateMultiConfig(Proc):
    """Generate a minimal cellranger multi config CSV (3' GEX only)
    from the downloaded FASTQ files.

    The output CSV is suitable for `cellranger multi --csv <outfile>`.
    Library rows use the sample-name prefix (text before the first '_S')
    and the directory that holds the FASTQs.
    """

    input = "fastqs:files"
    output = "outfile:file:multi_config.csv"
    input_data = lambda ch: [list(ch.iloc[:, 0])]
    requires = DownloadData
    envs = {"ref": config.ref.ref_cellranger_gex}
    lang = "python"
    script = r"""
from pathlib import Path
from panpath import LocalPath  # noqa: F401

fastqs = {{in.fastqs | each: as_path}}
ref    = {{envs.ref | repr}}
outfile = "{{out.outfile}}"

# All FASTQs land in the same directory
fastq_dir = str(Path(fastqs[0]).parent)

# Derive unique sample prefixes (text before the first "_S\d" token)
samples = []
seen = set()
for fq in fastqs:
    stem = Path(fq).name
    # e.g. "Sample_X_S1_L001_R1_001.fastq.gz" -> "Sample_X"
    prefix = stem.split("_S")[0]
    if prefix not in seen:
        seen.add(prefix)
        samples.append(prefix)

if not ref:
    raise ValueError(
        "Reference genome not configured. "
        "Set 'ref.ref_cellranger_gex' in your biopipen configuration "
        "(~/.biopipen.toml or ./.biopipen.toml)."
    )

with open(outfile, "w") as f:
    f.write("[gene-expression]\n")
    f.write(f"reference,{ref}\n")
    f.write("create-bam,false\n")
    f.write("\n")
    f.write("[libraries]\n")
    f.write("fastq_id,fastqs,feature_types\n")
    for sample in samples:
        f.write(f"{sample},{fastq_dir},Gene Expression\n")
"""

from pathlib import Path
class CellRangerMultiTest(CellRangerMulti):
    requires = CreateMultiConfig
    input_data = lambda ch: [(csv, None) for csv in ch.iloc[:, 0]]
    envs = {
        "ncores": 10,
        "cellranger": (
            "docker run -it --rm"
            f" -v {str(Path.cwd())}:{str(Path.cwd())}"
            " biopipen/cellranger:10.0.0"
        )
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
