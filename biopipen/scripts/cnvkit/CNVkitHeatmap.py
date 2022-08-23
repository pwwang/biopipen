from pathlib import Path

import cmdy
from diot import Diot

segfiles = {{in.segfiles | repr}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outdir = {{out.outdir | repr}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
by_bin = {{envs.by_bin | repr}}  # pyright: ignore
chromosome = {{ envs.chromosome | repr}}  # pyright: ignore
desaturate= {{ envs.desaturate | repr}}  # pyright: ignore
male_reference= {{ envs.male_reference | repr}}  # pyright: ignore
no_shift_xy= {{ envs.no_shift_xy | repr}}  # pyright: ignore
cases = {{envs.cases | repr}}  # pyright: ignore


def do_case(name, case):
    """Do heatmap with one case"""
    case = Diot(
        by_bin=by_bin,
        chromosome=chromosome,
        desaturate=desaturate,
        male_reference=male_reference,
        sample_sex=sample_sex or False,
        no_shift_xy=no_shift_xy,
        convert_args=convert_args,
    ) | case

    conv_args = case.pop("convert_args")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")
    cmdy.cnvkit.heatmap(
        **case,
        o=pdffile,
        _=segfiles,
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
