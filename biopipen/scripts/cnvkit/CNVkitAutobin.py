from __future__ import annotations

from pathlib import Path, PosixPath  # noqa: F401
from biopipen.utils.misc import run_command, dict_to_cli_args

bamfiles: list[Path] = {{in.bamfiles | each: as_path}}  # pyright: ignore # noqa
accfile: str | None = {{in.accfile | quote}}  # pyright: ignore
baitfile: str | None = {{in.baitfile | quote}}  # pyright: ignore
target_file: str | None = {{out.target_file | quote}}  # pyright: ignore
antitarget_file: str | None = {{out.antitarget_file | quote}}  # pyright: ignore
reffile: str = {{envs.ref | quote}}  # pyright: ignore
cnvkit: str = {{envs.cnvkit | quote}}  # pyright: ignore
method: str = {{envs.method | quote}}  # pyright: ignore
bp_per_bin: int = {{envs.bp_per_bin | repr}}  # pyright: ignore
target_max_size: int = {{envs.target_max_size | repr}}  # pyright: ignore
target_min_size: int = {{envs.target_min_size | repr}}  # pyright: ignore
antitarget_max_size: int = {{envs.antitarget_max_size | repr}}  # pyright: ignore
antitarget_min_size: int = {{envs.antitarget_min_size | repr}}  # pyright: ignore
annotate: str | None = {{envs.annotate | quote}}  # pyright: ignore
short_names: bool = {{envs.short_names | repr}}  # pyright: ignore

if baitfile == "None":
    baitfile = None
if accfile == "None":
    accfile = None


def main():

    args: dict = dict(
        f=Path(reffile).expanduser(),
        m=method,
        g=accfile,
        t=baitfile,
        b=bp_per_bin,
        target_max_size=target_max_size,
        target_min_size=target_min_size,
        antitarget_max_size=antitarget_max_size,
        antitarget_min_size=antitarget_min_size,
        annotate=False if annotate is None else Path(annotate).expanduser(),
        short_names=short_names,
        target_output_bed=target_file,
        antitarget_output_bed=antitarget_file,
        _=bamfiles,
    )
    args[""] = [cnvkit, "autobin"]

    run_command(dict_to_cli_args(args, dashify=True), fg=True)


if __name__ == "__main__":
    main()
