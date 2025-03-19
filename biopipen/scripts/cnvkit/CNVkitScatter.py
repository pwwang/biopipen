from pathlib import Path

from diot import Diot

from biopipen.utils.misc import run_command, dict_to_cli_args

cnrfile = {{in.cnrfile | quote}}  # pyright: ignore  # noqa
cnsfile = {{in.cnsfile | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
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
cases: dict | None = {{envs.cases | repr}}  # pyright: ignore


def do_case(name, case):
    """Do scatter with one case"""
    case = Diot(
        convert_args=convert_args,
        width=width,
        antitarget_marker=antitarget_marker,
        by_bin=by_bin,
        segment_color=False if segment_color is None else segment_color,
        trend=trend,
        y_max=False if y_max is None else y_max,
        y_min=False if y_min is None else y_min,
        min_variant_depth=min_variant_depth,
        zygosity_freq=zygosity_freq,
        title=False if title is None else title,
        chromosome=False if chromosome is None else chromosome,
        gene=False if gene is None else gene,
    ) | case

    conv_args = case.pop("convert_args")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")

    args: dict = dict(
        **case,
        s=cnsfile,
        o=pdffile,
        v=vcf,
        i=sample_id,
        n=normal_id,
        _=cnrfile,
    )
    args[""] = [cnvkit, "scatter"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

    conv_args: dict = dict(**conv_args, _=[pdffile, pngfile])
    conv_args[""] = [convert]
    run_command(
        dict_to_cli_args(conv_args, dashify=True, prefix="-"),
        fg=True,
    )


if __name__ == "__main__":
    if not cases:
        do_case("default", {})
    else:
        for name, case in cases.items():
            do_case(name, case)
