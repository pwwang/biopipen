from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args

excfiles = {{in.excfiles | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile = {{envs.ref | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
min_gap_size = {{envs.min_gap_size | quote}}  # pyright: ignore


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
