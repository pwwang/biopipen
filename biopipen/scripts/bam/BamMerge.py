from pathlib import Path
from biopipen.utils.misc import run_command

bamfiles = {{in.bamfiles | repr}}  # pyright: ignore
outfile = Path({{out.outfile | repr}})  # pyright: ignore
ncores = {{envs.ncores | int}}  # pyright: ignore
tool = {{envs.tool | quote}}  # pyright: ignore
samtools = {{envs.samtools | quote}}  # pyright: ignore
sambamba = {{envs.sambamba | quote}}  # pyright: ignore
should_sort = {{envs.sort | bool}}  # pyright: ignore
should_index = {{envs.index | bool}}  # pyright: ignore
merge_args = {{envs.merge_args | repr}}  # pyright: ignore
sort_args = {{envs.sort_args | repr}}  # pyright: ignore

if should_index and not should_sort:
    raise ValueError("Indexing requires sorting")


def use_samtools():
    """Use samtools to merge bam files"""
    print("Using samtools")
    ofile = (
        outfile
        if not should_sort
        else outfile.with_stem(f"{outfile.stem}.unsorted")
    )
    for key in ["-f", "-o", "-O", "--output-fmt", "-@", "--threads"]:
        if key in merge_args:
            raise ValueError(
                f"envs.merge_args cannot contain {key}, "
                "which is managed by the script"
            )
    cmd = [
        samtools,
        "merge",
        "-f",
        "-@",
        ncores,
        "-O",
        "BAM",
        "-o",
        ofile,
        *merge_args,
        *bamfiles,
    ]
    print("- Merging")
    run_command(cmd)

    if should_sort:
        print("- Sorting")
        for key in ["-o", "-@", "--threads"]:
            if key in sort_args:
                raise ValueError(
                    f"envs.sort_args cannot contain {key}, "
                    "which is managed by the script"
                )
        cmd = [
            samtools,
            "sort",
            "-@",
            ncores,
            *sort_args,
            "-o",
            outfile,
            ofile,
        ]
        run_command(cmd)

    if should_index:
        print("- Indexing")
        cmd = [samtools, "index", "-@", ncores, outfile]
        run_command(cmd)

    print("Done")


def use_sambamba():
    """Use sambamba to merge bam files"""
    print("Using sambamba")
    ofile = (
        outfile
        if not should_sort
        else outfile.with_stem(f"{outfile.stem}.unsorted")
    )
    for key in ["-t", "--nthreads"]:
        if key in merge_args:
            raise ValueError(
                f"envs.merge_args cannot contain {key}, "
                "which is managed by the script"
            )

    cmd = [sambamba, "merge", "-t", ncores, *merge_args, ofile, *bamfiles]
    print("- Merging")
    run_command(cmd)

    if should_sort:
        print("- Sorting")
        for key in ["-t", "--nthreads", "-o", "--out"]:
            if key in sort_args:
                raise ValueError(
                    f"envs.sort_args cannot contain {key}, "
                    "which is managed by the script"
                )

        cmd = [
            sambamba,
            "sort",
            "-t",
            ncores,
            *sort_args,
            "-o",
            outfile,
            ofile,
        ]
        run_command(cmd)

    if should_index:
        print("- Indexing")
        cmd = [sambamba, "index", "-t", ncores, outfile]
        run_command(cmd)

    print("Done")


if __name__ == "__main__":
    if tool == "samtools":
        use_samtools()
    elif tool == "sambamba":
        use_sambamba()
    else:
        raise ValueError(f"Unknown tool: {tool}")
