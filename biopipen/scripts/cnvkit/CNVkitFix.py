from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args

target_file = {{in.target_file | quote}}  # pyright: ignore  # noqa
antitarget_file = {{in.antitarget_file | quote}}  # pyright: ignore
reference_file = {{in.reference | quote}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
cluster = {{envs.cluster | repr}}  # pyright: ignore
no_gc = {{envs.no_gc | repr}}  # pyright: ignore
no_edge = {{envs.no_edge | repr}}  # pyright: ignore
no_rmask = {{envs.no_rmask | repr}}  # pyright: ignore

emptyfile = Path(outfile).parent.joinpath("empty")
emptyfile.touch()

def main():

    args: dict = dict(
        i=sample_id,
        o=outfile,
        c=cluster,
        no_gc=False,
        no_edge=False,
        no_rmask=False,
        _=[target_file, antitarget_file or emptyfile, reference_file],
    )
    args[""] = [cnvkit, "fix"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)


if __name__ == "__main__":
    main()
