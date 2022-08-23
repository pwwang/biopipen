from pathlib import Path

import cmdy

cnrfile = {{in.cnrfile | quote}}  # pyright: ignore
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
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

    cmdy.cnvkit.segment(
        o=outfile,
        d=Path(outfile).parent / "intermediate.rds",
        m=method,
        t=threshold,
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
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
