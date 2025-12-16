from pathlib import Path
from biopipen.utils.misc import logger

mquaddirs: list[str] = {{in.mquaddirs | each: str}}  # noqa: E999 # type: ignore
outdir: str = {{out.outdir | quote}}  # noqa: E999 # type: ignore


def merge_mtx_files(
    files: list[str],
    var_lists: list[list[str]],
    outfile: str,
    log_head: str,
) -> tuple[dict[str, int], list[str]]:

    tmpfile = outfile + ".tmp"
    ncells = {}
    all_vars = []
    with open(tmpfile, "w") as outfh:
        val_sum = 0
        prev_cells = 0
        for var_list, in_file in zip(var_lists, files):
            in_file = Path(in_file)
            sample = in_file.parent.stem
            logger.info(
                f"{log_head}: Merging {in_file.parent.stem} "
                f"with {len(var_list)} variants"
            )
            in_header = True
            with open(in_file, "r") as infh:
                for line in infh:
                    if line.startswith("%"):
                        continue
                    items = [int(x) for x in line.strip().split()]
                    if in_header:
                        if sample in ncells:
                            o = 1
                            while f"{sample}_{o}" in ncells:
                                o += 1

                            logger.warning(
                                f"{log_head}: Sample {sample} appears multiple times, "
                                f"rename to {sample}_{o}"
                            )
                            sample = f"{sample}_{o}"

                        ncells[sample] = items[1]
                        in_header = False
                        continue

                    val_sum += 1
                    variant = var_list[items[0] - 1]
                    if variant in all_vars:
                        var_idx = all_vars.index(variant) + 1
                    else:
                        all_vars.append(variant)
                        var_idx = len(all_vars)
                    outfh.write(f"{var_idx} {items[1] + prev_cells} {items[2]}\n")
            prev_cells += ncells[sample]

    with open(tmpfile, "r") as infh, open(outfile, "w") as outfh:
        outfh.write("%%MatrixMarket matrix coordinate integer general\n")
        outfh.write("%\n")
        outfh.write(f"{len(all_vars)} {sum(ncells.values())} {val_sum}\n")
        for line in infh:
            outfh.write(line)

    Path(tmpfile).unlink()
    return ncells, all_vars


var_lists = [
    [
        line.strip()
        for line in Path(f"{mquaddir}/passed_variant_names.txt")
        .read_text()
        .splitlines()
    ]
    for mquaddir in mquaddirs
]

out_ad_file = f"{outdir}/passed_ad.mtx"
merge_mtx_files(
    files=[f"{mquaddir}/passed_ad.mtx" for mquaddir in mquaddirs],
    var_lists=var_lists,
    outfile=out_ad_file,
    log_head="AD",
)

out_dp_file = f"{outdir}/passed_dp.mtx"
ncells, all_vars = merge_mtx_files(
    files=[f"{mquaddir}/passed_dp.mtx" for mquaddir in mquaddirs],
    var_lists=var_lists,
    outfile=out_dp_file,
    log_head="DP",
)

logger.info(f"Merging total {len(all_vars)} variants.")
out_pvn_file = f"{outdir}/passed_variant_names.txt"
Path(out_pvn_file).write_text("\n".join(all_vars) + "\n")

logger.info("Saving cell number info.")
out_cni_file = f"{outdir}/cell_number_info.txt"
with open(out_cni_file, "w") as outfh:
    outfh.write("Sample\tCellNumber\n")
    for sample, cellnum in ncells.items():
        outfh.write(f"{sample}\t{cellnum}\n")
