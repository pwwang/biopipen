from pathlib import Path
from biopipen.utils.misc import run_command, logger

# using:
# samtools view --subsample 0.1 --subsample-seed 1234 --threads 4 -b -o out.bam in.bam

bamfile = {{ in.bamfile | quote }} # pyright: ignore # noqa
bedfile = {{ in.bedfile | quote }} # pyright: ignore # noqa
outfile = Path({{ out.outfile | quote }}) # pyright: ignore
ncores = {{ envs.ncores | int }} # pyright: ignore
samtools = {{ envs.samtools | repr }} # pyright: ignore
tool = {{ envs.tool | repr }} # pyright: ignore
should_index = {{ envs.index | repr }} # pyright: ignore

if tool != "samtools":
    raise ValueError(
        f"Tool {tool} is not supported. "
        "Currently only samtools is supported."
    )

cmd = [
    samtools,
    "view",
    "--target-file",
    bedfile,
    "-b",
    "--threads",
    ncores,
    "-o",
    outfile,
    bamfile
]
run_command(cmd, fg=True)

if should_index:
    logger.info("Indexing the output bam file.")
    cmd = [samtools, "index", "-@", ncores, outfile]
    run_command(cmd, fg=True)
