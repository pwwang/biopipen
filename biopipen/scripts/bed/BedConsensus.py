import sys
from math import ceil
from pathlib import Path

import cmdy

bedfiles = {{in.bedfiles | repr}}  # pyright: ignore
outfile = Path({{out.outbed | repr}})  # pyright: ignore
bedtools_path = {{envs.bedtools | repr}}  # pyright: ignore
cutoff = {{envs.cutoff | repr}}  # pyright: ignore
distance = {{envs.distance | repr}}  # pyright: ignore
chrsize = {{envs.chrsize | repr}}  # pyright: ignore
bedfiles = [Path(bedfile) for bedfile in bedfiles]
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
    cmdy.bedtools.sort(i=concatfile, _exe=bedtools_path).r() > sortedfile

    return sortedfile


def genomecov():
    """Calculate genome coverage."""
    _log("- Calculating genome coverage")
    genomecovfile = outfile.parent / "_genomecov.bed"
    filteredfile = outfile.parent / "_filtered.bed"
    cmdy.bedtools.genomecov(
        i=concat_bedfiles(),
        g=chrsize,
        bg=True,
        _exe=bedtools_path,
        _prefix="-",
    ).r() > genomecovfile
    cmdy.awk(f'$4 >= {cutoff}', genomecovfile).r() > filteredfile
    return filteredfile


def merge_if_needed(awfile):
    _log("- Merging results if needed using distance:", distance)
    if distance == 0:
        awfile.rename(outfile)
        return

    cmdy.bedtools.merge(
        i=awfile,
        d=distance,
        c=4,
        o="collapse",
        _exe=bedtools_path,
    ).r() > outfile


def main():
    filteredfile = genomecov()
    merge_if_needed(filteredfile)


if __name__ == "__main__":
    main()
