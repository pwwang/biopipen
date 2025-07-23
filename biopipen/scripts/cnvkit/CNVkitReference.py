from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args

covfiles = {{in.covfiles | default: [] | each: str | repr}}  # pyright: ignore  # noqa
target_file = {{in.target_file | quote: quote_none=False}}  # pyright: ignore
antitarget_file = {{in.antitarget_file | quote: quote_none=False}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile: str = {{envs.ref | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
cluster = {{envs.cluster | repr}}  # pyright: ignore
min_cluster_size = {{envs.min_cluster_size | repr}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}} # pyright: ignore
no_gc = {{envs.no_gc | repr}}  # pyright: ignore
no_edge = {{envs.no_edge | repr}}  # pyright: ignore
no_rmask = {{envs.no_rmask | repr}}  # pyright: ignore


def main():

    args = dict(
        f=Path(reffile).expanduser(),
        o=outfile,
        c=cluster,
        min_cluster_size=min_cluster_size,
        x=sample_sex,
        y=male_reference,
        t=False if not target_file or covfiles else target_file,
        a=False if not antitarget_file or covfiles else antitarget_file,
        no_gc=False,
        no_edge=False,
        no_rmask=False,
        _=covfiles or [],
    )
    args[""] = [cnvkit, "reference"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)


if __name__ == "__main__":
    main()
