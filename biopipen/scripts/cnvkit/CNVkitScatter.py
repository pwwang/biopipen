from pathlib import Path

import cmdy
from diot import Diot

cnrfile = {{in.cnrfile | quote}}  # pyright: ignore
cnsfile = {{in.cnsfile | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
chromosome = {{envs.chromosome | repr}}  # pyright: ignore
gene = {{envs.gene | repr}}  # pyright: ignore
width = {{envs.width | repr}}  # pyright: ignore
antitarget_marker = {{envs.antitarget_marker | repr}}  # pyright: ignore
by_bin = {{envs.by_bin | repr}}  # pyright: ignore
segment_color = {{envs.segment_color | repr}}  # pyright: ignore
trend = {{envs.trend | repr}}  # pyright: ignore
y_max = {{envs.y_max | repr}}  # pyright: ignore
y_min = {{envs.y_min | repr}}  # pyright: ignore
min_variant_depth = {{envs.min_variant_depth | repr}}  # pyright: ignore
zygosity_freq = {{envs.zygosity_freq | repr}}  # pyright: ignore
title = {{envs.title | repr}}  # pyright: ignore
cases = {{envs.cases | repr}}  # pyright: ignore


def do_case(name, case):
    """Do scatter with one case"""
    case = Diot(
        convert_args=convert_args,
        width=width,
        antitarget_marker=antitarget_marker,
        by_bin=by_bin,
        segment_color=segment_color,
        trend=trend,
        y_max=y_max,
        y_min=y_min,
        min_variant_depth=min_variant_depth,
        zygosity_freq=zygosity_freq,
        title=title,
        chromosome=chromosome,
        gene=gene,
    ) | case

    conv_args = case.pop("convert_args")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")

    cmdy.cnvkit.scatter(
        **case,
        s=cnsfile,
        o=pdffile,
        v=vcf or False,
        i=sample_id or False,
        n=normal_id or False,
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
