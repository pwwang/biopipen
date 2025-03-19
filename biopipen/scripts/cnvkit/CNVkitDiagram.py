from pathlib import Path

from diot import Diot
from biopipen.utils.misc import run_command, dict_to_cli_args


cnrfile = {{in.cnrfile | quote}}  # pyright: ignore  # noqa
cnsfile = {{in.cnsfile | quote}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
threshold = {{envs.threshold | repr}}  # pyright: ignore
min_probes = {{envs.min_probes | repr}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}}  # pyright: ignore
no_shift_xy = {{envs.no_shift_xy | repr}}  # pyright: ignore
title = {{envs.title | repr}}  # pyright: ignore
cases: dict | None = {{envs.cases | repr}}  # pyright: ignore


def do_case(name, case):
    """Do scatter with one case"""
    case = Diot(
        convert_args=convert_args,
        threshold=threshold,
        min_probes=min_probes,
        male_reference=male_reference,
        sample_sex=sample_sex or False,
        no_shift_xy=no_shift_xy,
        title=False if title is None else title,
    ) | case

    conv_args = case.pop("convert_args")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")

    args: dict = dict(
        **case,
        s=cnsfile,
        o=pdffile,
        _=cnrfile,
    )
    args[""] = [cnvkit, "diagram"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

    conv_args: dict = dict(**conv_args, _=[pdffile, pngfile])
    conv_args[""] = [convert]
    run_command(
        dict_to_cli_args(conv_args, prefix="-", dashify=True),
        fg=True,
    )


if __name__ == "__main__":
    if not cases:
        do_case("default", {})
    else:
        for name, case in cases.items():
            do_case(name, case)
