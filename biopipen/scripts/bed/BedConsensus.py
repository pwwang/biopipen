from pathlib import Path

import cmdy
import pandas

from datar import f
from datar.base import as_double, mean, paste
from datar.dplyr import (
    mutate,
    bind_rows,
    group_by,
    summarise,
    filter_,
    select,
    ungroup,
)

bedfiles = {{in.bedfiles | repr}}  # pyright: ignore
outfile = Path({{out.outbed | repr}})  # pyright: ignore
bedtools_path = {{envs.bedtools | repr}}  # pyright: ignore
binsize = {{envs.binsize | repr}}  # pyright: ignore
ignore_scores = {{envs.ignore_scores | repr}}  # pyright: ignore
cutoff = {{envs.cutoff | repr}}  # pyright: ignore
distance = {{envs.distance | repr}}  # pyright: ignore
chrsize = {{envs.chrsize | repr}}  # pyright: ignore
bedfiles = [Path(bedfile) for bedfile in bedfiles]


def read_bed(bedfile, bedidx):
    """Read BED file."""
    ofile = outfile.parent / f"_{bedfile.stem}.bed"
    df = pandas.read_csv(bedfile, sep="\t", header=None)
    header = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    df.columns = header[:len(df.columns)]
    if "score" in df.columns and bedidx not in ignore_scores:
        # average to single base
        df = df >> mutate(score=f.score / as_double(f.end - f.start))
    else:
        df = df >> mutate(score=f.end - f.start)

    df.to_csv(ofile, sep="\t", index=False, header=False)
    return ofile


def makewindows():
    """Make windows using binsize."""
    cmdy.bedtools.makewindows(
        g=chrsize,
        w=binsize,
        _exe=bedtools_path,
    ).r() > (outfile.parent / "windows.bed")

    return outfile.parent / "windows.bed"


def get_weights(windowfile, bedfile, bedidx):
    """Get weights."""
    ofile = outfile.parent / f"_{bedfile.stem}_binned.bed"
    owfile = outfile.parent / f"_{bedfile.stem}_weighted_bin.bed"

    bedfile = read_bed(bedfile, bedidx)

    cmdy.bedtools.intersect(
        a=windowfile,
        b=bedfile,
        wo=True,
        _prefix="-",
        _exe=bedtools_path,
    ).r() > ofile

    df = pandas.read_csv(ofile, sep="\t", header=None)
    header = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ]
    df.columns = header[:len(df.columns) - 1] + ["overlap"]
    df = df >> mutate(
        weight=as_double(f.score) *
        (f.overlap / binsize) /
        (f.end2 - f.start2)
    )
    df.to_csv(owfile, sep="\t", index=False, header=True)

    return owfile


def avg_weights_and_filter(owfiles):
    ofile = outfile.parent / "_avg_weights_filtered.bed"
    df = None
    for owfile in owfiles:
        tmp = pandas.read_csv(owfile, sep="\t", header=0)
        df = df >> bind_rows(tmp)

    df = df >> group_by(f.chrom1, f.start1, f.end1) >> summarise(
        chrom=f.chrom1,
        start=f.start1,
        end=f.end1,
        name=paste(f.name, collapse=":"),
        score=mean(f.weight),
        strand="+",
    ) >> filter_(
        f.score >= cutoff
    ) >> ungroup() >> select(
        ~f.chrom1, ~f.start1, ~f.end1,
    )

    df.to_csv(ofile, sep="\t", index=False, header=False)
    return ofile, len(df.columns)


def merge_if_needed(awfile, ncols):
    if distance == 0:
        cmdy.cp(awfile, outfile)
    else:
        cmdy.bedtools.merge(
            i=awfile,
            d=distance,
            c=f"{ncols-2},{ncols-1}",
            o="distinct,mean",
            _exe=bedtools_path,
        ).r() > outfile


def main():
    binfile = makewindows()
    owfiles = [
        get_weights(binfile, bedfile, i)
        for i, bedfile in enumerate(bedfiles)
    ]
    awfile, ncols = avg_weights_and_filter(owfiles)
    merge_if_needed(awfile, ncols)


if __name__ == "__main__":
    main()
