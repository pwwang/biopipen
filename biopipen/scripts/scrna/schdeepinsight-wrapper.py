from argparse import ArgumentParser
import os
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

from SCHdeepinsight import immune  # noqa: E402

classifier = immune(args.outdir)
results = classifier.run_pipeline(
    input_file=args.input,
    ref_file=args.ref,
    batch_size=args.batch_size
)
results.to_csv(args.output, sep="\t", index=False)
