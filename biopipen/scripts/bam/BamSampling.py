from pathlib import Path
from biopipen.utils.misc import run_command, logger

# using:
# samtools view --subsample 0.1 --subsample-seed 1234 --threads 4 -b -o out.bam in.bam

bamfile = {{ in.bamfile | quote }} # pyright: ignore # noqa
outfile = Path({{ out.outfile | quote }}) # pyright: ignore
ncores = {{ envs.ncores | int }} # pyright: ignore
samtools = {{ envs.samtools | repr }} # pyright: ignore
tool = {{ envs.tool | repr }} # pyright: ignore
fraction: float = {{ envs.fraction | repr }} # pyright: ignore
seed = {{ envs.seed | int }} # pyright: ignore
should_index = {{ envs.index | repr }} # pyright: ignore
should_sort = {{ envs.sort | repr }} # pyright: ignore
sort_args = {{ envs.sort_args | repr }} # pyright: ignore

if should_index and not should_sort:
    raise ValueError("Indexing requires sorting")

if fraction is None:
    raise ValueError("'envs.fraction' must be provided.")

if tool != "samtools":
    raise ValueError(
        f"Tool {tool} is not supported. "
        "Currently only samtools is supported."
    )

if fraction > 1:
    # calculate the fraction based on the number of reads
    logger.info("Converting fraction > 1 to a fraction of reads.")
    cmd = [
        samtools,
        "view",
        "--threads",
        ncores,
        "-c",
        bamfile
    ]
    nreads = run_command(cmd, stdout="return").strip()  # type: ignore
    fraction = fraction / float(int(nreads))

ofile = (
    outfile
    if not should_sort
    else outfile.with_stem(f"{outfile.stem}.unsorted")
)

cmd = [
    samtools,
    "view",
    "--subsample",
    fraction,
    "--subsample-seed",
    seed,
    "--threads",
    ncores,
    "-b",
    "-o",
    ofile,
    bamfile
]
run_command(cmd, fg=True)

if should_sort:
    logger.info("Sorting the output bam file.")
    for key in ["-o", "-@", "--threads"]:
        if key in sort_args:
            raise ValueError(
                f"envs.sort_args cannot contain {key}, "
                "which is managed by the script"
            )

    cmd = [
        samtools,
        "sort",
        "-@",
        ncores,
        *sort_args,
        "-o",
        outfile,
        ofile
    ]
    run_command(cmd, fg=True)

if should_index:
    logger.info("Indexing the output bam file.")
    cmd = [samtools, "index", "-@", ncores, outfile]
    run_command(cmd, fg=True)
