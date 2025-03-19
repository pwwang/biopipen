from pathlib import Path

import pandas
from biopipen.utils.misc import run_command, dict_to_cli_args

metafile: str = {{in.cnsfile | quote}}  # pyright: ignore  # noqa
outdir: str = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
method = {{envs.method | quote}}  # pyright: ignore
segment_method = {{envs.segment_method | quote}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}}  # pyright: ignore
count_reads = {{envs.count_reads | repr}}  # pyright: ignore
drop_low_coverage = {{envs.drop_low_coverage | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
rscript = {{envs.rscript | quote}}  # pyright: ignore
ref: str = {{envs.ref | quote}}  # pyright: ignore
targets = {{envs.targets | repr}}  # pyright: ignore
antitargets = {{envs.antitargets | repr}}  # pyright: ignore
annotate = {{envs.annotate | repr}}  # pyright: ignore
short_names = {{envs.short_names | repr}}  # pyright: ignore
target_avg_size = {{envs.target_avg_size | repr}}  # pyright: ignore
access = {{envs.access | repr}}  # pyright: ignore
access_min_gap_size = {{envs.access_min_gap_size | repr}}  # pyright: ignore
access_excludes = {{envs.access_excludes | repr}}  # pyright: ignore
antitarget_avg_size = {{envs.antitarget_avg_size | repr}}  # pyright: ignore
antitarget_min_size = {{envs.antitarget_min_size | repr}}  # pyright: ignore
cluster = {{envs.cluster | repr}}  # pyright: ignore
reference = {{envs.reference | repr}}  # pyright: ignore
scatter = {{envs.scatter | repr}}  # pyright: ignore
diagram = {{envs.diagram | repr}}  # pyright: ignore
type_tumor = {{envs.type_tumor | repr}}  # pyright: ignore
type_normal = {{envs.type_normal | repr}}  # pyright: ignore
type_col: str = {{envs.type_col | quote}}  # pyright: ignore


def gen_access():
    if access:
        return access

    accessfile = Path(outdir) / "access.bed"
    args: dict = dict(
        exclude=access_excludes or False,
        s=access_min_gap_size or False,
        o=accessfile,
        _=Path(ref).expanduser(),
    )
    args[""] = [cnvkit, "access"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

    return accessfile


def parse_metafile():
    """Parse the metafile to get the tumor (case) and normal (control)
    bam files."""
    meta = pandas.read_csv(metafile, sep="\t", header=0)
    if "BamFile" not in meta.columns:
        raise ValueError("BamFile column not found in metafile")
    if type_col not in meta.columns:
        raise ValueError(f"{type_col} column not found in metafile")

    tumor_bams = meta.loc[meta[type_col] == type_tumor, "BamFile"].tolist()
    normal_bams = (
        meta
        .loc[meta[type_col] == type_normal, "BamFile"]
        .tolist()
    )

    if len(tumor_bams) == 0:
        raise ValueError(f"No {type_tumor} samples found in metafile")

    return tumor_bams, normal_bams


def main():
    accessfile = gen_access()
    tumor_bams, normal_bams = parse_metafile()

    args = dict(
        method=method,
        segment_method=segment_method,
        y=male_reference,
        c=count_reads,
        drop_low_coverage=drop_low_coverage,
        p=ncores,
        rscript_path=rscript,
        n=normal_bams or False,
        f=ref,
        t=targets or False,
        a=antitargets or False,
        annotate=annotate,
        short_names=short_names,
        target_avg_size=target_avg_size,
        access=accessfile,
        antitarget_avg_size=antitarget_avg_size,
        antitarget_min_size=antitarget_min_size,
        cluster=cluster,
        reference=reference,
        scatter=scatter,
        diagram=diagram,
        d=outdir,
        _=tumor_bams,
    )
    args[""] = [cnvkit, "batch"]
    run_command(dict_to_cli_args(args, dashify=True), fg=True)


if __name__ == "__main__":
    main()
