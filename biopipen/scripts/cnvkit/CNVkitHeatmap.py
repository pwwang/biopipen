from pathlib import Path

from diot import Diot  # type: ignore[import]

from biopipen.utils.misc import run_command, dict_to_cli_args

segfiles = {{in.segfiles | default: [] | each: str | repr}}  # pyright: ignore # noqa  # noqa
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
convert = {{envs.convert | quote}}  # pyright: ignore
convert_args = {{envs.convert_args | repr}}  # pyright: ignore
by_bin = {{envs.by_bin | repr}}  # pyright: ignore
chromosome = {{ envs.chromosome | repr}}  # pyright: ignore
desaturate= {{ envs.desaturate | repr}}  # pyright: ignore
male_reference= {{ envs.male_reference | repr}}  # pyright: ignore
no_shift_xy= {{ envs.no_shift_xy | repr}}  # pyright: ignore
order = {{envs.order | repr}}  # pyright: ignore
cases: dict | None = {{envs.cases | repr}}  # pyright: ignore


def parse_order(files, orderfile):
    if not orderfile:
        return files

    def _match(path, sample):
        p = Path(path)
        return sample in (p.name, p.stem)

    out = []
    files = set(files)
    orders = Path(orderfile).read_text().splitlines()
    for od in orders:
        matched = None
        for fpath in files:
            if _match(fpath, od):
                matched = fpath
                break
        if matched:
            out.append(matched)
            files.remove(matched)

    return out


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
        order=order,
    ) | case

    conv_args = case.pop("convert_args")
    the_order = case.pop("order")
    pdffile = Path(outdir).joinpath(f"{name}.heatmap.pdf")
    pngfile = Path(outdir).joinpath(f"{name}.heatmap.png")

    sfiles = parse_order(segfiles, the_order)
    args = dict(
        **case,
        o=pdffile,
        _=sfiles,
    )
    args[""] = [cnvkit, "heatmap"]
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
