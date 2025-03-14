import sys
from math import ceil
from pathlib import Path, PosixPath  # noqa: F401

from biopipen.utils.misc import run_command

bedfiles: list[Path] = {{in.bedfiles | each: as_path}}  # pyright: ignore # noqa
outfile = Path({{out.outbed | quote}})  # pyright: ignore
bedtools_path = {{envs.bedtools | repr}}  # pyright: ignore
cutoff: float = {{envs.cutoff | repr}}  # pyright: ignore
distance = {{envs.distance | repr}}  # pyright: ignore
chrsize = {{envs.chrsize | repr}}  # pyright: ignore
# bedfiles = [Path(bedfile) for bedfile in bedfiles]
# In case there are duplicated stems
stems = [f"{bedfile.stem}__{i}" for i, bedfile in enumerate(bedfiles)]

if cutoff < 1:
    cutoff = ceil(len(bedfiles) * cutoff)


def _log(*msg):
    print(*[str(m) for m in msg], file=sys.stderr)


def concat_bedfiles():
    """Concatenate and merge bedfiles."""
    concatfile = outfile.parent / "_concatenated.bed"
    sortedfile = outfile.parent / "_sorted.bed"
    concatfile.write_text("")

    _log("- Concatenating BED files")
    with open(concatfile, "a") as fout:
        for bedfile in bedfiles:
            fout.write(bedfile.read_text())

    _log("- Sorting the concatenated BED file")
    run_command(
        [bedtools_path, "sort", "-i", concatfile, "-faidx", chrsize],
        stdout=sortedfile,
    )

    return sortedfile


def genomecov():
    """Calculate genome coverage."""
    _log("- Calculating genome coverage")
    genomecovfile = outfile.parent / "_genomecov.bed"
    filteredfile = outfile.parent / "_filtered.bed"
    run_command(
        [
            bedtools_path,
            "genomecov",
            "-i",
            concat_bedfiles(),
            "-g",
            chrsize,
            "-bg",
        ],
        stdout=genomecovfile,
    )
    run_command(
        ["awk", f'$4 >= {cutoff}', genomecovfile],
        stdout=filteredfile,
    )
    return filteredfile


def merge_if_needed(awfile):
    _log("- Merging results if needed using distance:", distance)
    if distance == 0:
        awfile.rename(outfile)
        return

    run_command(
        [
            bedtools_path,
            "merge",
            "-i",
            awfile,
            "-d",
            distance,
            "-c",
            "4",
            "-o",
            "collapse",
        ],
        stdout=outfile,
    )


def main():
    filteredfile = genomecov()
    merge_if_needed(filteredfile)


if __name__ == "__main__":
    main()
