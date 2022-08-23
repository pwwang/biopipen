from pathlib import Path
import cmdy

cnsfile = {{in.cnsfile | quote}}  # pyright: ignore
cnrfile = {{in.cnrfile | quote}}  # pyright: ignore
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
center = {{envs.center | repr}}  # pyright: ignore
center_at = {{envs.center_at | repr}}  # pyright: ignore
filter = {{envs.filter | repr}}  # pyright: ignore
method = {{envs.method | quote}}  # pyright: ignore
thresholds = {{envs.thresholds | repr}}  # pyright: ignore
ploidy = {{envs.ploidy | repr}}  # pyright: ignore
purity = {{envs.purity | repr}}  # pyright: ignore
drop_low_coverage = {{envs.drop_low_coverage | repr}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}}  # pyright: ignore
min_variant_depth = {{envs.min_variant_depth | repr}}  # pyright: ignore
zygosity_freq = {{envs.zygosity_freq | repr}}  # pyright: ignore


def main():
    outfile = (
        Path(outdir).joinpath(Path(cnsfile).name).with_suffix(".call.cns")
    )
    vcffile = Path(outdir).joinpath(Path(cnsfile).name).with_suffix(".vcf")
    bedfile = Path(outdir).joinpath(Path(cnsfile).name).with_suffix(".bed")

    cmdy.cnvkit.call(
        # thresholds could be "-1.1...", starting `-` make it unregconizable
        f"-t={thresholds}" if thresholds else "",
        center=center,
        center_at=center_at,
        filter=filter,
        m=method,
        ploidy=ploidy,
        purity=purity,
        drop_low_coverage=drop_low_coverage,
        male_reference=male_reference,
        sample_sex=sample_sex or False,
        o=outfile,
        v=vcf or False,
        i=sample_id or False,
        n=normal_id or False,
        min_variant_depth=min_variant_depth,
        zygosity_freq=zygosity_freq,
        _=cnsfile,
        _exe=cnvkit,
    ).fg()

    # conver to vcf
    cmdy.cnvkit.export(
        "vcf",
        cnr=cnrfile,
        i=sample_id or False,
        ploidy=ploidy,
        sample_sex=sample_sex or False,
        male_reference=male_reference,
        o=vcffile,
        _=cnsfile,
        _exe=cnvkit,
    ).fg()

    # convert to bed
    cmdy.cnvkit.export(
        "bed",
        i=sample_id or False,
        ploidy=ploidy,
        sample_sex=sample_sex or False,
        male_reference=male_reference,
        label_genes=True,
        show="ploidy",  # all, variant, ploidy
        o=bedfile,
        _=cnsfile,
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
