from hashlib import md5
from pathlib import Path
from biopipen.utils.misc import run_command, dict_to_cli_args

infile: str = {{ in.bamfile | quote }} # pyright: ignore # noqa
outfile = Path({{ out.outfile | quote }}) # pyright: ignore
args: dict = {{ envs | dict | repr }} # pyright: ignore
ncores = args.pop("ncores")
tool = args.pop("tool")
samtools = args.pop("samtools")
sambamba = args.pop("sambamba")
tmpdir = args.pop("tmpdir")
byname = args.pop("byname")
should_index = args.pop("index")
sig = md5(infile.encode()).hexdigest()
tmpdir = Path(tmpdir).joinpath(
    f"biopipen_BamSort_{{job.index}}_{sig}_{Path(infile).name}"
)
tmpdir.mkdir(parents=True, exist_ok=True)
tmpdir = str(tmpdir)


def use_samtools():
    """Use samtools to sort/index bam file.

    Usage: samtools sort [options...] [in.bam]
    Options:
    -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
    -u         Output uncompressed data (equivalent to -l 0)
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
    -M         Use minimiser for clustering unaligned/unplaced reads
    -K INT     Kmer size to use for minimiser [20]
    -n         Sort by read name (not compatible with samtools index command)
    -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
    -o FILE    Write final output to FILE rather than standard output
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam
        --no-PG
                Do not add a PG line
        --template-coordinate
                Sort by template-coordinate
        --input-fmt-option OPT[=VAL]
                Specify a single input file format option in the form
                of OPTION or OPTION=VALUE
    -O, --output-fmt FORMAT[,OPT[=VAL]]...
                Specify output format (SAM, BAM, CRAM)
        --output-fmt-option OPT[=VAL]
                Specify a single output file format option in the form
                of OPTION or OPTION=VALUE
        --reference FILE
                Reference sequence FASTA FILE [null]
    -@, --threads INT
                Number of additional threads to use [0]
        --write-index
                Automatically index the output files [off]
        --verbosity INT
                Set level of verbosity
    """  # noqa
    sargs = args.copy()
    sargs["n"] = byname
    sargs["T"] = f"{tmpdir}/tmp"
    sargs["threads"] = ncores

    if should_index:
        sargs["write-index"] = True
        # https://github.com/samtools/samtools/issues/1196
        sargs["o"] = f"{outfile}##idx##{outfile}.bai"
    else:
        sargs["o"] = outfile

    n_outfmt = sum(["O" in sargs, "output-fmt" in sargs])
    if n_outfmt > 1:
        raise ValueError(
            "envs.args cannot contain both 'O' and 'output-fmt'"
        )
    if n_outfmt == 0:
        sargs["O"] = "BAM"

    cmd = [
        samtools,
        "sort",
        *dict_to_cli_args(sargs),
        infile,
    ]
    run_command(cmd)


def use_sambamba():
    """Use sambamba to sort/index bam file.

    sambamba 0.8.2
    by Artem Tarasov and Pjotr Prins (C) 2012-2021
        LDC 1.28.1 / DMD v2.098.1 / LLVM12.0.0 / bootstrap LDC - the LLVM D compiler (1.28.1)

    Usage: sambamba-sort [options] <input.bam>

    Options: -m, --memory-limit=LIMIT
                approximate total memory limit for all threads (by default 2GB)
            --tmpdir=TMPDIR
                directory for storing intermediate files; default is system directory for temporary files
            -o, --out=OUTPUTFILE
                output file name; if not provided, the result is written to a file with .sorted.bam extension
            -n, --sort-by-name
                sort by read name instead of coordinate (lexicographical order)
            --sort-picard
                sort by query name like in picard
            -N, --natural-sort
                sort by read name instead of coordinate (so-called 'natural' sort as in samtools)
            -M, --match-mates
                pull mates of the same alignment together when sorting by read name
            -l, --compression-level=COMPRESSION_LEVEL
                level of compression for sorted BAM, from 0 to 9
            -u, --uncompressed-chunks
                write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space
            -p, --show-progress
                show progressbar in STDERR
            -t, --nthreads=NTHREADS
                use specified number of threads
            -F, --filter=FILTER
                keep only reads that satisfy FILTER
    """  # noqa
    sargs = args.copy()
    sargs["nthreads"] = ncores
    sargs["n"] = byname
    sargs["tmpdir"] = tmpdir
    sargs["o"] = outfile
    cmd = [
        sambamba,
        "sort",
        *dict_to_cli_args(sargs, sep="="),
        infile,
    ]
    run_command(cmd)


if __name__ == "__main__":
    if tool == "samtools":
        use_samtools()
    elif tool == "sambamba":
        use_sambamba()
    else:
        raise ValueError(f"Unknown tool: {tool}")
