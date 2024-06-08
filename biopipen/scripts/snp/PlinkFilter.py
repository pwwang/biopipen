"""Script for snp.PlinkFilter"""

from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

indir = {{in.indir | repr}}  # pyright: ignore # noqa: #999
samples_file_keep = {{in.samples_file_keep | repr}}  # pyright: ignore
variants_file_keep = {{in.variants_file_keep | repr}}  # pyright: ignore
samples_file_remove = {{in.samples_file_remove | repr}}  # pyright: ignore
variants_file_remove = {{in.variants_file_remove | repr}}  # pyright: ignore
outdir = {{out.outdir | repr}}  # pyright: ignore

plink = {{envs.plink | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
samples_keep = {{envs.samples_keep | repr}}  # pyright: ignore
variants_keep = {{envs.variants_keep | repr}}  # pyright: ignore
samples_remove = {{envs.samples_remove | repr}}  # pyright: ignore
variants_remove = {{envs.variants_remove | repr}}  # pyright: ignore
e_samples_file_keep = {{envs.samples_file_keep | repr}}  # pyright: ignore
e_variants_file_keep = {{envs.variants_file_keep | repr}}  # pyright: ignore
e_samples_file_remove = {{envs.samples_file_remove | repr}}  # pyright: ignore
e_variants_file_remove = {{envs.variants_file_remove | repr}}  # pyright: ignore
chr = {{envs.chr | repr}}  # pyright: ignore
not_chr = {{envs.not_chr | repr}}  # pyright: ignore
autosome = {{envs.autosome | repr}}  # pyright: ignore
autosome_xy = {{envs.autosome_xy | repr}}  # pyright: ignore
snps_only = {{envs.snps_only | repr}}  # pyright: ignore

samples_file_keep = samples_file_keep or e_samples_file_keep
if not samples_file_keep and samples_keep:
    samples_file_keep = Path(outdir) / "samples_keep.txt"
    if isinstance(samples_keep, str):
        samples_keep = [s.strip() for s in samples_keep.split(",")]

    with open(samples_file_keep, "w") as fh:
        fh.writelines(samples_keep)

variants_file_keep = variants_file_keep or e_variants_file_keep
if not variants_file_keep and variants_keep:
    variants_file_keep = Path(outdir) / "variants_keep.txt"
    if isinstance(variants_keep, str):
        variants_keep = [v.strip() for v in variants_keep.split(",")]

    with open(variants_file_keep, "w") as fh:
        fh.writelines(variants_keep)

samples_file_remove = samples_file_remove or e_samples_file_remove
if not samples_file_remove and samples_remove:
    samples_file_remove = Path(outdir) / "samples_remove.txt"
    if isinstance(samples_remove, str):
        samples_remove = [s.strip() for s in samples_remove.split(",")]

    with open(samples_file_remove, "w") as fh:
        fh.writelines(samples_remove)

variants_file_remove = variants_file_remove or e_variants_file_remove
if not variants_file_remove and variants_remove:
    variants_file_remove = Path(outdir) / "variants_remove.txt"
    if isinstance(variants_remove, str):
        variants_remove = [v.strip() for v in variants_remove.split(",")]

    with open(variants_file_remove, "w") as fh:
        fh.writelines(variants_remove)

bedfile = list(Path(indir).glob("*.bed"))
if len(bedfile) == 0:
    raise FileNotFoundError(f"No .bed file found in `in.indir`")
elif len(bedfile) > 1:
    logger.warning(f"Multiple .bed files found in `in.indir`, using the first one.")

bedfile = bedfile[0]
input = bedfile.with_suffix("")
output = Path(outdir) / bedfile.stem

args = {
    "": [plink],
    "bfile": input,
    "out": output,
    "threads": ncores,
    "keep_allele_order": True,
    "make-bed": True,
}

if samples_file_keep:
    args["keep"] = samples_file_keep
if variants_file_keep:
    args["extract"] = variants_file_keep
if samples_file_remove:
    args["remove"] = samples_file_remove
if variants_file_remove:
    args["exclude"] = variants_file_remove
if chr:
    args["chr"] = chr
if not_chr:
    args["not_chr"] = not_chr
if autosome:
    args["autosome"] = True
if autosome_xy:
    args["autosome"] = True
if snps_only:
    args["snps_only"] = snps_only

run_command(dict_to_cli_args(args, dashify=True), fg=True)
