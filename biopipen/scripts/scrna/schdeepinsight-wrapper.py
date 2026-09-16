"""Run scHDeepInsight on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-schdeepinsight.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
scHDeepInsight is gated by default: the `SCHdeepinsight` package (0.3.5) is not
installed in this environment, and the tool additionally needs a reference file
(`--ref`, the immune-cell `reference.rds` of the repo) plus the pretrained model
checkpoint from <https://github.com/shangruJia/scHDeepInsight>. To un-gate it,
install `pip install SCHdeepinsight` (and
`pip install git+https://github.com/alok-ai-lab/pyDeepInsight.git`) into the
python running this wrapper and download the checkpoint from the repo.
The import below reports all of this when the package is missing, instead of a
bare ModuleNotFoundError.
"""
from argparse import ArgumentParser
import os
import sys
import pandas as pd

parser = ArgumentParser(description="Run scHDeepInsight")
parser.add_argument(
    "-i", "--input", required=True, help="Input H5AD file (AnnData)"
)
parser.add_argument("-r", "--ref", required=True, help="Reference RDS file")
parser.add_argument(
    "-o", "--output", required=True, help="Output TSV file for results"
)
parser.add_argument(
    "-d", "--outdir", required=True,
    help="Output directory for intermediate files"
)
parser.add_argument(
    "-b", "--batch-size", type=int, default=128,
    help="Batch size for CNN prediction"
)
parser.add_argument(
    "--rhome", required=True, help="R home directory (for rpy2)"
)

args = parser.parse_args()

# Set R_HOME for rpy2 before importing scHDeepInsight
os.environ["R_HOME"] = args.rhome

try:
    from SCHdeepinsight import immune  # noqa: E402
except ImportError as exc:
    sys.exit(
        f"Cannot import scHDeepInsight ({exc}). The `SCHdeepinsight` package "
        "(0.3.5) is not installed in the python running this wrapper, and the "
        "tool also needs a reference file (`--ref`) and the pretrained model "
        "checkpoint from https://github.com/shangruJia/scHDeepInsight, so it is "
        "gated by default.\n"
        "To un-gate it, into that python (`envs.schdeepinsight.python`):\n"
        "  pip install SCHdeepinsight\n"
        "  pip install git+https://github.com/alok-ai-lab/pyDeepInsight.git\n"
        "and download the checkpoint from the repo."
    )

classifier = immune(args.outdir)
results = classifier.run_pipeline(
    input_file=args.input,
    ref_file=args.ref,
    batch_size=args.batch_size
)
results.to_csv(args.output, sep="\t", index=False)
