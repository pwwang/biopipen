from pathlib import Path, PosixPath  # noqa: F401
from biopipen.utils.misc import run_command, dict_to_cli_args

excfiles: list[Path] = {{in.excfiles | each: as_path}}  # pyright: ignore # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore
reffile: str = {{envs.ref | quote}}  # pyright: ignore
cnvkit: str = {{envs.cnvkit | quote}}  # pyright: ignore
min_gap_size: str = {{envs.min_gap_size | quote}}  # pyright: ignore


def main():
    other_args = {
        "": [cnvkit, "access"],
        "s": min_gap_size,
        "o": outfile,
        "_": Path(reffile).expanduser(),
    }
    if excfiles:
        other_args["exclude"] = excfiles

    run_command(dict_to_cli_args(other_args, dashify=True), fg=True)


if __name__ == "__main__":
    main()
