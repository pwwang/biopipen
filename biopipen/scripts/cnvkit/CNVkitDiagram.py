from pathlib import Path

import cmdy
from diot import Diot


cnrfile = {{in.cnrfile | quote}}  # pyright: ignore
cnsfile = {{in.cnsfile | quote}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
threshold = {{envs.threshold | repr}}  # pyright: ignore
min_probes = {{envs.min_probes | repr}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}}  # pyright: ignore
no_shift_xy = {{envs.no_shift_xy | repr}}  # pyright: ignore
title = {{envs.title | repr}}  # pyright: ignore
cases = {{envs.cases | repr}}  # pyright: ignore


def do_case(name, case):
    """Do scatter with one case"""
    case = Diot(
        convert_args=convert_args,
        threshold=threshold,
        min_probes=min_probes,
        male_reference=male_reference,
        sample_sex=sample_sex or False,
        no_shift_xy=no_shift_xy,
        title=title,
    ) | case

    conv_args = case.pop("convert_args")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")

    cmdy.cnvkit.diagram(
        **case,
        s=cnsfile,
        o=pdffile,
        _=cnrfile,
        _exe=cnvkit,
    ).fg()

    cmdy.convert(
        **conv_args,
        _=[pdffile, pngfile],
        _prefix="-",
        _exe=convert,
    ).fg()


if __name__ == "__main__":
    if not cases:
        do_case("default", {})
    else:
        for name, case in cases.items():
            do_case(name, case)
