import sys
from multiprocessing import Pool
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
ncores = {{envs.ncores | repr}}  # pyright: ignore
bedfiles = [Path(bedfile) for bedfile in bedfiles]
# In case there are duplicated stems
stems = [f"{bedfile.stem}__{i}" for i, bedfile in enumerate(bedfiles)]


def _log(*msg):
    print(*[str(m) for m in msg], file=sys.stderr)


def read_bed(bedfile, bedidx):
    """Read BED file."""
    _log("- Reading BED file:", bedfile)
    ofile = outfile.parent / f"_{stems[bedidx]}.bed"
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
        ofile = bedfile
    else:
        df = df >> mutate(score=f.end - f.start)
        df.to_csv(ofile, sep="\t", index=False, header=False)

    return ofile


def makewindows():
    """Make windows using binsize."""
    _log("- Making windows with binsize:", binsize)
    cmdy.bedtools.makewindows(
        g=chrsize,
        w=binsize,
        _exe=bedtools_path,
    ).r() > (outfile.parent / "windows.bed")

    return outfile.parent / "windows.bed"


def get_weights(windowfile, bedfile, bedidx):
    """Get weights."""
    _log("- Getting bin weights for:", bedfile)
    ofile = outfile.parent / f"_{stems[bedidx]}_binned.bed"
    owfile = outfile.parent / f"_{stems[bedidx]}_weighted_bin.bed"

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
    _log("- Averaging bin weights")
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


def avg_weights_and_filter_fast(owfiles):
    """A faster version of avg_weights_and_filter."""
    _log("- Averaging bin weights")
    ofile = outfile.parent / "_avg_weights_filtered.bed"
    mfile = outfile.parent / "_avg_weights_merged.bed"
    sfile = outfile.parent / "_avg_weights_merged_sorted.bed"
    # merge the files
    header = None
    with open(mfile, "w") as fm:
        for owfile in owfiles:
            with open(owfile, "r") as fow:
                for i, line in enumerate(fow):
                    if i > 0:
                        fm.write(line)
                    elif header is None:
                        header = line

    cmdy.sort(
        _prefix="-",
        _=mfile,
        **{"k1,1": True, "k2,2n": True}
    ).r() > sfile

    # loop and aggregate
    names = []
    scores = []
    key = None
    with open(sfile, "r") as fs, open(ofile, "w") as fo:
        for line in fs:
            parts = line.rstrip("\n").split("\t")
            my_key = "\t".join(parts[:3])
            if my_key != key and key is not None:
                name = ":".join(names)
                score = sum(scores) / len(scores)
                if score >= cutoff:
                    fo.write(f"{key}\t{name}\t{score}\t+\n")
                names = [parts[6]]
                scores = [float(parts[-1])]
            else:
                names.append(parts[6])
                scores.append(float(parts[-1]))
            key = my_key
        name = ":".join(names)
        score = sum(scores) / len(scores)
        fo.write(f"{key}\t{name}\t{score}\t+\n")

    return ofile, 6


def merge_if_needed(awfile, ncols):
    _log("- Merging results if needed using distance:", distance)
    if distance == 0:
        cmdy.cp(awfile, outfile)
    else:
        awfile_sorted = outfile.parent / f"{awfile.stem}_sorted{awfile.suffix}"
        cmdy.sort(
            _prefix="-",
            _=awfile,
            **{"k1,1": True, "k2,2n": True}
        ).r() > awfile_sorted
        cmdy.bedtools.merge(
            i=awfile_sorted,
            d=distance,
            c=f"{ncols-2},{ncols-1}",
            o="distinct,mean",
            _exe=bedtools_path,
        ).r() > outfile


def main():
    binfile = makewindows()

    with Pool(ncores) as pool:
        owfiles = pool.starmap(
            get_weights,
            [(binfile, bedfile, i) for i, bedfile in enumerate(bedfiles)],
        )

    awfile, ncols = avg_weights_and_filter_fast(owfiles)
    merge_if_needed(awfile, ncols)


if __name__ == "__main__":
    main()
