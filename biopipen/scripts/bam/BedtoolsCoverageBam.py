from __future__ import annotations

from biopipen.utils.misc import run_command, dict_to_cli_args

bamfile: str = {{ in.bamfile | str | quote }}  # pyright: ignore  # noqa
bedfile: str = {{ in.bedfile | str | quote }}  # pyright: ignore  # noqa
outfile: str = {{ out.outfile | str | quote }}  # pyright: ignore

envs: dict = {{ envs | dict | repr }}  # pyright: ignore
bedtools: str = envs.pop("bedtools", "bedtools")

envs["a"] = bedfile
envs["b"] = bamfile

cmd = [
    bedtools,
    "coverage",
    *dict_to_cli_args(envs, prefix="-", dup_key=False)
]
run_command(cmd, stdout=outfile)
