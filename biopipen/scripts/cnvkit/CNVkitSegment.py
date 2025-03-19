from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args

cnrfile = {{in.cnrfile | quote}}  # pyright: ignore  # noqa
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
method = {{envs.method | quote}}  # pyright: ignore
threshold = {{envs.threshold | repr}}  # pyright: ignore
drop_low_coverage = {{envs.drop_low_coverage | repr}}  # pyright: ignore
drop_outliers = {{envs.drop_outliers | repr}}  # pyright: ignore
rscript = {{envs.rscript | quote}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
smooth_cbs = {{envs.smooth_cbs | repr}}  # pyright: ignore
min_variant_depth = {{envs.min_variant_depth | repr}}  # pyright: ignore
zygosity_freq = {{envs.zygosity_freq | repr}}  # pyright: ignore


def main():

    args: dict = dict(
        o=outfile,
        d=Path(outfile).parent / "intermediate.rds",
        m=method,
        t=False if threshold is None else threshold,
        drop_low_coverage=drop_low_coverage,
        drop_outliers=drop_outliers,
        rscript_path=rscript,
        p=ncores,
        smooth_cbs=smooth_cbs,
        v=vcf or False,
        i=sample_id or False,
        n=normal_id or False,
        min_variant_depth=min_variant_depth,
        zygosity_freq=zygosity_freq,
        _=cnrfile,
    )
    args[""] = [cnvkit, "segment"]
    cmd_args = dict_to_cli_args(args, dashify=True)
    run_command(cmd_args, fg=True)


if __name__ == "__main__":
    main()
