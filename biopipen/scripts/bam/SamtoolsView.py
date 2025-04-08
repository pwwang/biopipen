from pathlib import PosixPath  # type: ignore # noqa
from biopipen.utils.misc import run_command, dict_to_cli_args
from biopipen.utils.reference import bam_index

bamfile: str = {{ in.bamfile | quote }}  # pyright: ignore # noqa
outfile: str = {{ out.outfile | quote }}  # pyright: ignore # noqa
envs: dict = {{envs | attr: "to_dict" | call}}  # pyright: ignore  # noqa
ncores = envs.pop("ncores")
samtools = envs.pop("samtools")
should_index = envs.pop("index")


def run_samtools(infile):
    cmd = [
        samtools,
        "view",
        "-b",
        "--threads",
        str(ncores),
        "-o",
        outfile,
    ] + dict_to_cli_args(envs, dashify=True) + [infile]

    run_command(cmd, fg=True)
    if should_index:
        bam_index(outfile, tool="samtools", samtools=samtools, ncores=ncores)

    return outfile


if __name__ == "__main__":
    infile = bam_index(bamfile, tool="samtools", samtools=samtools, ncores=ncores)
    run_samtools(infile)
