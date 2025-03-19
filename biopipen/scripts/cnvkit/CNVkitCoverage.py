from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args

bamfile: str = {{in.bamfile | quote}}  # pyright: ignore  # noqa
target_file = {{in.target_file | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile: str = {{envs.reffile | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
count = {{envs.count | repr}}  # pyright: ignore
min_mapq = {{envs.min_mapq | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore


def main():

    args = dict(
        f=Path(reffile).expanduser(),
        c=count,
        q=min_mapq,
        p=ncores,
        o=outfile,
        _=[bamfile, target_file],
    )
    args[""] = [cnvkit, "coverage"]
    run_command(dict_to_cli_args(args), fg=True)


if __name__ == "__main__":
    main()
