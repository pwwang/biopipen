from __future__ import annotations

from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

indir: str = {{in.indir | quote}}  # pyright: ignore # noqa: #999
samples_file = {{in.samples_file | quote}}  # pyright: ignore
variants_file = {{in.variants_file | quote}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore

plink = {{envs.plink | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
samples: list[str] | str = {{envs.samples | repr}}  # pyright: ignore
variants: list[str] | str = {{envs.variants | repr}}  # pyright: ignore
e_samples_file = {{envs.samples_file | repr}}  # pyright: ignore
e_variants_file = {{envs.variants_file | repr}}  # pyright: ignore
keep = {{envs.keep | repr}}  # pyright: ignore
vfile_type = {{envs.vfile_type | repr}}  # pyright: ignore
chr = {{envs.chr | repr}}  # pyright: ignore
not_chr = {{envs.not_chr | repr}}  # pyright: ignore
autosome = {{envs.autosome | repr}}  # pyright: ignore
autosome_xy = {{envs.autosome_xy | repr}}  # pyright: ignore
snps_only = {{envs.snps_only | repr}}  # pyright: ignore

samples_file = samples_file or e_samples_file
if not samples_file and samples:
    samples_file = Path(outdir) / "_samples.txt"
    if isinstance(samples, str):
        samples = [s.strip() for s in samples.split(",")]

    with open(samples_file, "w") as fh:
        fh.writelines(
            [
                line.replace("/", "\t") + "\n"
                if "/" in line
                else line + "\t" + line + "\n"
                for line in samples
            ]
        )

variants_file = variants_file or e_variants_file
if not variants_file and variants:
    if vfile_type != "id":
        logger.warning(
            "envs.vfile_type should be 'id' if only envs.variants is provided."
        )
        vfile_type = "id"

    variants_file = Path(outdir) / "_variants.txt"
    if isinstance(variants, str):
        variants = [v.strip() for v in variants.split(",")]

    with open(variants_file, "w") as fh:
        fh.writelines([line + "\n" for line in variants])

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
    "make-bed": True,
}

if keep:
    if samples_file:
        args["keep"] = samples_file
    if variants_file:
        args["extract"] = (
            variants_file if vfile_type == "id" else [vfile_type, variants_file]
        )
else:
    if samples_file:
        args["remove"] = samples_file
    if variants_file:
        args["exclude"] = (
            variants_file if vfile_type == "id" else [vfile_type, variants_file]
        )

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

run_command(dict_to_cli_args(args, dashify=True, dup_key=False), fg=True)
