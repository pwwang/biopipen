from pathlib import Path
from biopipen.utils.misc import run_command

cnsfile: str = {{in.cnsfile | quote}}  # pyright: ignore  # noqa
cnrfile: str = {{in.cnrfile | quote}}  # pyright: ignore
vcf = {{in.vcf | repr}}  # pyright: ignore
sample_id = {{in.sample_id | repr}}  # pyright: ignore
normal_id = {{in.normal_id | repr}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
purity = {{in.purity | repr}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
center = {{envs.center | repr}}  # pyright: ignore
center_at = {{envs.center_at | repr}}  # pyright: ignore
filter = {{envs.filter | repr}}  # pyright: ignore
method = {{envs.method | quote}}  # pyright: ignore
thresholds = {{envs.thresholds | repr}}  # pyright: ignore
ploidy = {{envs.ploidy | repr}}  # pyright: ignore
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

    cmd = [cnvkit, "call"]
    cmd.extend(["--center", center])
    cmd.extend(["--method", method])
    cmd.extend(["--ploidy", ploidy])
    cmd.extend(["--output", outfile])
    if filter is not None:
        cmd.extend(["--filter", filter])
    if center_at is not None:
        cmd.extend(["--center-at", center_at])
    if thresholds and method == "threshold":
        cmd.append(f"--thresholds={thresholds}")
    if purity is not None:
        cmd.extend(["--purity", purity])
    if drop_low_coverage:
        cmd.append("--drop-low-coverage")
    if male_reference:
        cmd.append("--male-reference")
    if sample_sex:
        cmd.extend(["--sample-sex", sample_sex])
    if min_variant_depth is not None:
        cmd.extend(["--min-variant-depth", min_variant_depth])
    if zygosity_freq is True:
        cmd.append("--zygosity-freq")
    elif zygosity_freq is not None:
        cmd.extend(["--zygosity-freq", zygosity_freq])
    if vcf:
        cmd.extend(["--vcf", vcf])
    if sample_id:
        cmd.extend(["--sample-id", sample_id])
    if normal_id:
        cmd.extend(["--normal-id", normal_id])
    cmd.append(cnsfile)

    run_command(cmd, fg=True)

    # conver to vcf
    # If this errored, see:
    # https://github.com/etal/cnvkit/pull/818
    cmd = [cnvkit, "export", "vcf"]
    cmd.extend(["--ploidy", ploidy])
    cmd.extend(["--output", vcffile])
    if sample_id:
        cmd.extend(["--sample-id", sample_id])
    if sample_sex:
        cmd.extend(["--sample-sex", sample_sex])
    if male_reference:
        cmd.append("--male-reference")
    cmd.append(cnsfile)

    run_command(cmd, fg=True)

    # convert to bed
    cmd = [cnvkit, "export", "bed"]
    cmd.extend(["--ploidy", ploidy])
    cmd.extend(["--output", bedfile])
    cmd.append("--label-genes")
    cmd.extend(["--show", "ploidy"])  # all, variant, ploidy
    if sample_id:
        cmd.extend(["--sample-id", sample_id])
    if sample_sex:
        cmd.extend(["--sample-sex", sample_sex])
    if male_reference:
        cmd.append("--male-reference")
    cmd.append(cnsfile)

    run_command(cmd, fg=True)


if __name__ == "__main__":
    main()
