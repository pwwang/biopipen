#!/usr/bin/env python
"""Generate tiny synthetic test data for the splitpipe tests

The test data is generated locally because `split-pipe` (Parse Biosciences
Pipeline) is required by the tests but its example data and reference genomes
must not be downloaded. The generator only reads the per-round barcode tables
(`bc_data_*.csv`) of an installed/bundled split-pipe; those tables are NEVER
copied into this repo, they are read at generation time from `--barcodes-dir`
(or from the `barcodes/` directory of the installed `splitpipe` package).

Barcode design (verified against split-pipe v1.9.0):
    Each R2 read is `UMI(10) + bc3 + L12 + bc2 + L23 + bc1` (86 bp, the
    chemistry-v2 amplicon from `splitpipe.kits.get_chem_amp_seq`); split-pipe
    slices exactly those offsets. Each barcode must come from the set of its
    OWN round — for `--chemistry v2 --kit WT` split-pipe's kit table assigns
    bc1=n99_v5, bc2=v1, bc3=v1 (see `ROUND_FILES`) — because split-pipe
    checks every field against its round's set and then counts perfect
    combinations; barcodes from the wrong round leave that count empty and the
    run dies in `bcutils.get_perfect_bc_counts` ("argmin of an empty
    sequence").

    bc1 is baked in per cell, and the two sublibraries use DISJOINT sets of
    bc1 barcodes, so no cell barcode is shared between the sublibraries when
    they are combined. bc2/bc3 vary per cell (each cell is one barcode
    combination) and the per-read UMI makes every R2 unique.

Outputs (all in `--outdir`, defaulting to this directory):
    ref.fa            synthetic reference, 4 contigs (`chr1`..`chr4`),
                      ~75 kb each by default
    ref.gtf           gene annotation, 3 genes per contig (12 genes), 2-4
                      exons per gene, spliced transcripts of ~1.5-3 kb, so
                      reads sampled from spliced transcripts span exon-exon
                      junctions
    mkref.par         parameter file for `split-pipe --mode mkref` with a
                      small `genomeSAindexNbases` for STAR, passed through
                      the single `mkref_star_args` parameter as a
                      comma-separated list, since STAR requires
                      `genomeSAindexNbases <= log2(genome_length) / 2 - 1`,
                      which is far below the STAR default of 14 for a genome
                      this small
    <s1,s2>_R1.fastq.gz   cDNA reads (70-90 bp) sampled from the spliced
                      transcripts, quality fixed to "I"
    <s1,s2>_R2.fastq.gz   the four 8-mer barcodes of the read, quality fixed
                      to "I"

Header convention (required by split-pipe):
    The read names follow the Illumina convention
    `@<read-name> 1:N:0:<index>` in R1 and `@<read-name> 2:N:0:<index>` in
    R2, with the same `<read-name>` for the two mates of a pair. split-pipe's
    `fastq_head_get_name_num` (`spc_setup.py:608`) splits each header on " "
    (or "/") and requires the second part to start with "1" (R1) or "2" (R2).
    Headers without that mate field make `--mode all` fail while setting up
    the fastqs ("Failed to find <name> <num> separater").

Verified against a real run of split-pipe v1.9.0 (`--mode all --chemistry v2
--kit WT`):
    - the header convention above, with the shared read name and the `1`/`2`
      mate field: this is what lets `--mode all` set up both fastqs
    - `--chemistry v2` and `--kit WT`: the run log shows "Set bc_amp_seq for
      chemistry 2" and "Set kit to WT, given"

Usage:

    python tests/test_splitpipe/data/make_test_data.py \\
        --barcodes-dir /path/to/splitpipe/barcodes

The script asserts that no generated file reaches 1 MB.
"""
from __future__ import annotations

import argparse
import gzip
import csv
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).parent
BASES = "ACGT"
MAX_FILE_SIZE = 1_000_000  # 1 MB

CONTIGS = 4
GENES_PER_CONTIG = 3
READ_LEN_RANGE = (70, 90)
R2_FIELDS = 4
BC_LEN = 8


def parse_args() -> argparse.Namespace:
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--barcodes-dir",
        default=None,
        help=(
            "The directory with the barcode tables (bc_dict_<chemistry>.json "
            "and bc_data_<chemistry>.csv) of a Parse Biosciences Pipeline "
            "installation. Defaults to the `barcodes/` directory of the "
            "installed `splitpipe` package."
        ),
    )
    parser.add_argument(
        "--chemistry",
        choices=["v1", "v2"],
        default="v2",
        help="The chemistry version selecting bc_dict_<chemistry>.json",
    )
    parser.add_argument(
        "--outdir",
        default=str(HERE),
        help="The output directory of the generated data",
    )
    parser.add_argument("--seed", type=int, default=42, help="The random seed")
    parser.add_argument(
        "--cells",
        type=int,
        default=120,
        help="The number of cells per sublibrary",
    )
    parser.add_argument(
        "--reads-per-cell",
        type=int,
        default=40,
        help="The number of read pairs per cell",
    )
    parser.add_argument(
        "--contig-size",
        type=int,
        default=75000,
        help="The length of each contig",
    )
    return parser.parse_args()


def resolve_barcodes_dir(cli_dir: str | None) -> Path:
    """Resolve the barcode table directory

    Either from `--barcodes-dir` or from the installed `splitpipe` package.
    """
    if cli_dir:
        barcode_dir = Path(cli_dir)
    else:
        try:
            import splitpipe
        except ImportError:
            sys.exit(
                "Could not import the `splitpipe` package to locate its "
                "barcode tables. Pass --barcodes-dir pointing at the "
                "`barcodes/` directory of a Parse Biosciences Pipeline "
                "installation, e.g. `--barcodes-dir "
                "/path/to/splitpipe/barcodes`."
            )
        barcode_dir = Path(splitpipe.__file__).parent / "barcodes"
        if not barcode_dir.is_dir():
            sys.exit(
                f"The `barcodes/` directory was not found in the installed "
                f"`splitpipe` package: {barcode_dir}. Pass --barcodes-dir "
                f"pointing at the `barcodes/` directory of a Parse "
                f"Biosciences Pipeline installation."
            )

    if not barcode_dir.is_dir():
        sys.exit(f"The barcode directory does not exist: {barcode_dir}")
    return barcode_dir


# The barcode rounds split-pipe uses for `--chemistry v2 --kit WT`, from its own
# kit table (`kits.KIT_CHEM_DTAB`, printed by `kits.kit_bc_set_list("WT", "v2")`):
# bc1 = n99_v5, bc2 = v1, bc3 = v1. Note the file *named* `bc_data_v2.csv` belongs
# to chemistry 1 (`--chemistry v1`), not to `--chemistry v2`. The fixture is pinned
# to the v2/WT pair, so the round files are pinned here too — re-derive them with
# `kits.kit_bc_set_list("WT", <chem>)` if either changes.
ROUND_FILES = {
    "bc1": "bc_data_n99_v5.csv",
    "bc2": "bc_data_v1.csv",
    "bc3": "bc_data_v1.csv",
}


def load_round_barcodes(barcode_dir: Path) -> dict[str, list[str]]:
    """Load the valid 8-mer barcode set of each R2 barcode round

    Every emitted R2 must draw its three barcodes from these sets: split-pipe
    checks each field against the set of its round and then counts perfect
    combinations. If none exist, `bcutils.get_perfect_bc_counts` dies with
    `ValueError: attempt to get argmin of an empty sequence`.
    """
    sets: dict[str, list[str]] = {}
    for round_name, fname in ROUND_FILES.items():
        table = barcode_dir / fname
        if not table.is_file():
            sys.exit(f"The barcode table does not exist: {table}")
        with table.open() as fh:
            seqs = sorted({row["sequence"] for row in csv.DictReader(fh)})
        if not seqs or any(len(bc) != BC_LEN for bc in seqs):
            sys.exit(f"Unexpected barcodes in {table}: {seqs[:5]} ...")
        sets[round_name] = seqs
    return sets


def build_reference(
    rng: random.Random,
    contig_size: int,
) -> tuple[dict[str, str], list[dict]]:
    """Build the synthetic contigs and genes

    Returns the contig sequences (name -> sequence) and the genes, each with
    the exons as 0-based half-open intervals on its contig and the spliced
    transcript sequence.
    """
    contigs: dict[str, str] = {}
    genes: list[dict] = []
    for ci in range(CONTIGS):
        ctg = f"chr{ci + 1}"
        seq = "".join(rng.choice(BASES) for _ in range(contig_size))
        contigs[ctg] = seq
        slot = contig_size // GENES_PER_CONTIG
        for gi in range(GENES_PER_CONTIG):
            n_exons = rng.randint(2, 4)
            tx_len = rng.randint(1500, 3000)
            exon_lens = [tx_len // n_exons] * n_exons
            for i in range(tx_len - exon_lens[0] * n_exons):
                exon_lens[i] += 1

            pos = gi * slot + rng.randint(200, 400)
            exons: list[tuple[int, int]] = []
            for i, exon_len in enumerate(exon_lens):
                if i:
                    pos += rng.randint(300, 1500)  # intron
                exons.append((pos, pos + exon_len))
                pos += exon_len
            if pos > (gi + 1) * slot:
                sys.exit(
                    f"Gene {gi + 1} of {ctg} does not fit into its slot; "
                    "use a larger --contig-size."
                )

            genes.append(
                {
                    "name": f"GENE{ci * GENES_PER_CONTIG + gi + 1}",
                    "gid": f"gene{ci * GENES_PER_CONTIG + gi + 1}",
                    "tid": f"tx{ci * GENES_PER_CONTIG + gi + 1}",
                    "contig": ctg,
                    "exons": exons,
                    "seq": "".join(seq[start:end] for start, end in exons),
                }
            )
    return contigs, genes


def write_fasta(contigs: dict[str, str], path: Path, width: int = 60) -> None:
    """Write the contigs to a fasta file"""
    with path.open("w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(f"{seq[i : i + width]}\n")


def write_gtf(genes: list[dict], path: Path) -> None:
    """Write the genes to a GTF file

    Coordinates are converted from 0-based half-open intervals to 1-based
    inclusive ones. Features are written in coordinate order, as expected by
    STAR.
    """
    with path.open("w") as fh:
        for gene in genes:
            gid, tid, gname = gene["gid"], gene["tid"], gene["name"]
            attrs = (
                f'gene_id "{gid}"; transcript_id "{tid}"; '
                f'gene_name "{gname}"; gene_biotype "protein_coding";'
            )
            tx_start = gene["exons"][0][0] + 1
            tx_end = gene["exons"][-1][1]
            fh.write(
                f"{gene['contig']}\tsynthetic\tgene\t{tx_start}\t{tx_end}\t"
                f".\t+\t.\t{attrs}\n"
            )
            fh.write(
                f"{gene['contig']}\tsynthetic\ttranscript\t{tx_start}\t"
                f"{tx_end}\t.\t+\t.\t{attrs}\n"
            )
            for start, end in gene["exons"]:
                fh.write(
                    f"{gene['contig']}\tsynthetic\texon\t{start + 1}\t{end}\t"
                    f".\t+\t.\t{attrs}\n"
                )


def write_parfile(genome_len: int, path: Path) -> int:
    """Write the `split-pipe --mode mkref` parameter file

    STAR requires `genomeSAindexNbases <= log2(genome_length) / 2 - 1`, which
    is ~8 for the tiny genome generated here, far below the STAR default of
    14. STAR options are passed through the single `mkref_star_args`
    parameter as a comma-separated list, because split-pipe lowercases
    parfile keys and does not recognize the camelCase
    `mkref_star_genomeSAindexNbases` of Parse's shipped example
    (`CRISPR-focal_spipe.par`).
    """
    nbases = max(1, min(14, int(math.log2(genome_len) / 2 - 1)))
    with path.open("w") as fh:
        fh.write("# Parameter file for `split-pipe --mode mkref`\n")
        fh.write(
            f"# The synthetic genome is {genome_len} bp long. STAR's "
            f"suffix-array index needs\n"
        )
        fh.write(
            f"# genomeSAindexNbases <= log2({genome_len}) / 2 - 1 = "
            f"{math.log2(genome_len) / 2 - 1:.2f}, hence {nbases}.\n"
        )
        fh.write(
            "# STAR options go through the single `mkref_star_args` parameter "
            "as a\n"
        )
        fh.write(
            "# comma-separated list: split-pipe lowercases parfile keys, so "
            "the camelCase\n"
        )
        fh.write(
            "# `mkref_star_genomeSAindexNbases` of Parse's shipped example is "
            "not recognized.\n"
        )
        fh.write(
            f"mkref_star_args  --genomeSAindexNbases,{nbases},"
            f"--genomeSAsparseD,2,--limitGenomeGenerateRAM,4000000000\n"
        )
    return nbases


def write_count_parfile(path: Path, cells: int) -> None:
    """Write the `split-pipe --mode all` parameter file

    Each fixture sublibrary has only 4,800 reads, far below split-pipe's
    100,000-read minimum estimate (`fq_size_guess_min`), which would abort
    before any alignment.
    """
    with path.open("w") as fh:
        fh.write(
            "# Parameter file for `split-pipe --mode all` on the tiny "
            "fixture data.\n"
        )
        fh.write(
            "# Each sublibrary has only 4,800 reads, far below split-pipe's "
            "100,000-read\n"
        )
        fh.write(
            "# minimum estimate (`fq_size_guess_min`), which would abort "
            "before alignment.\n"
        )
        fh.write("fq_size_guess_min  1000\n")
        fh.write(
            "# Tiny-data overrides, each verified against split-pipe v1.9.0.\n"
            "# `fastq_samp_slice` keeps the sample real: `check_fastq_samp_slice`\n"
            "# only rewrites the slice while guess*guess_fac > its end, and with\n"
            "# end <= 0.75 X 4,800 the slice survives — the default end lets the\n"
            "# start be clamped to read 100,000, which makes the sample empty and\n"
            "# crashes the chemistry scoring / InQC. `inqc_min_reads` is InQC's\n"
            "# read-count minimum (default 100,000) and the two `*_score_skip`\n"
            "# flags drop the chemistry/kit scoring, which needs real-scale data.\n"
        )
        fh.write("fastq_samp_slice  0 3000\n")
        fh.write("inqc_min_reads  1000\n")
        fh.write("chem_score_skip  True\n")
        fh.write("kit_score_skip  True\n")
        fh.write(
            "# Remaining floors, all from split-pipe's own parameter list:\n"
            "# multiprocess chunk splitting indexes one record per jump and\n"
            "# needs more than the tiny file provides; the alignment wants\n"
            "# 50,000 reads; DGE's TSCP/cell cuts assume real-scale data, and\n"
            "# the 'given' hints keep it from float()-ing an unset estimate.\n"
        )
        fh.write("pre_num_multiproc  1\n")
        fh.write("align_min_start_reads  1000\n")
        fh.write("dge_min_start_tscp  100\n")
        fh.write("dge_tscp_cut_min  3\n")
        fh.write("dge_tscp_cut_given  3\n")
        fh.write("dge_cell_min  2\n")
        fh.write("dge_cell_hard_min  2\n")
        fh.write(f"dge_cell_given  {cells}\n")
        fh.write("pre_perfect_bc_cell_min  2\n")


# Chemistry v2 R2 layout, from split-pipe's own v2 amplicon
# (`splitpipe.kits.get_chem_amp_seq("v2")`):
#   NNNNNNNNNN 33333333 GTGGCCGATGTTTCGCATCGGCGTACGACT 22222222
#   ATCCACGTGCTTGAGACTGTGG 11111111
# i.e. UMI(10) + barcode3 + L12 + barcode2 + L23 + barcode1. split-pipe
# slices exactly these offsets (`bcutils.get_amp_part_bounds`), so the
# barcodes must sit here and the linkers must be present. Re-derive if
# Parse changes the chemistry definition.
UMI_LEN = 10
LINKER12 = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
LINKER23 = "ATCCACGTGCTTGAGACTGTGG"


def make_reads(
    rng: random.Random,
    sublib: str,
    f1_barcodes: list[str],
    bc_sets: dict[str, list[str]],
    genes: list[dict],
    cells: int,
    reads_per_cell: int,
    seen_r2: set[str],
    seen_cells: set[tuple[str, str, str]],
) -> list[tuple[str, str, str]]:
    """Make the reads of one sublibrary

    Returns a list of (read name, R1 sequence, R2 sequence). R2 follows the
    chemistry-v2 layout verified against split-pipe: a 10-bp UMI, then the
    cell's three barcodes separated by the L12/L23 linkers. Each barcode comes
    from the set of its OWN round (`bc_sets`: bc1 from n99_v5, bc2/bc3 from v1
    for kit WT) — the combination identifies the cell; bc1 alone is disjoint
    between sublibraries; the UMI makes every read unique.
    """
    reads: list[tuple[str, str, str]] = []
    for ci in range(cells):
        bc1 = f1_barcodes[ci % len(f1_barcodes)]
        while True:
            bc2 = rng.choice(bc_sets["bc2"])
            bc3 = rng.choice(bc_sets["bc3"])
            if (bc1, bc2, bc3) not in seen_cells:
                seen_cells.add((bc1, bc2, bc3))
                break

        for ri in range(reads_per_cell):
            gene = rng.choice(genes)
            tx = gene["seq"]
            read_len = rng.randint(*READ_LEN_RANGE)
            start = rng.randrange(len(tx) - read_len + 1)
            r1 = tx[start : start + read_len]
            while True:
                umi = "".join(rng.choices("ACGT", k=UMI_LEN))
                r2 = umi + bc3 + LINKER12 + bc2 + LINKER23 + bc1
                if r2 not in seen_r2:
                    seen_r2.add(r2)
                    break
            reads.append((f"{sublib}_cell{ci}_read{ri}", r1, r2))
    return reads


def write_fastqs(
    reads: list[tuple[str, str, str]],
    fq1: Path,
    fq2: Path,
) -> None:
    """Write the read pairs to two gzipped fastq files

    The headers follow the Illumina convention `@<name> 1:N:0:<index>` (R1)
    and `@<name> 2:N:0:<index>` (R2), with the same `<name>` for the two
    mates, as required by split-pipe's `fastq_head_get_name_num`
    (`spc_setup.py:608`).

    `mtime=0` keeps the compressed bytes identical between runs.
    """
    with (
        gzip.GzipFile(fq1, "wb", mtime=0) as fh1,
        gzip.GzipFile(fq2, "wb", mtime=0) as fh2,
    ):
        for name, r1, r2 in reads:
            fh1.write(f"@{name} 1:N:0:ACGTACGT\n{r1}\n+\n{'I' * len(r1)}\n".encode())
            fh2.write(f"@{name} 2:N:0:ACGTACGT\n{r2}\n+\n{'I' * len(r2)}\n".encode())


def main() -> None:
    """Generate the test data"""
    args = parse_args()
    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    barcode_dir = resolve_barcodes_dir(args.barcodes_dir)
    bc_sets = load_round_barcodes(barcode_dir)
    print(
        f"# Loaded {sum(len(v) for v in bc_sets.values())} {BC_LEN}-mer barcodes for chemistry "
        f"{args.chemistry} from {barcode_dir}"
    )

    contigs, genes = build_reference(rng, args.contig_size)
    genome_len = sum(len(seq) for seq in contigs.values())

    outputs = [
        outdir / "ref.fa",
        outdir / "ref.gtf",
        outdir / "mkref.par",
        outdir / "s1_R1.fastq.gz",
        outdir / "s1_R2.fastq.gz",
        outdir / "s2_R1.fastq.gz",
        outdir / "s2_R2.fastq.gz",
    ]
    write_fasta(contigs, outputs[0])
    write_gtf(genes, outputs[1])
    nbases = write_parfile(genome_len, outputs[2])
    write_count_parfile(outputs[2].parent / "count.par", args.cells)
    print(
        f"# Reference: {len(contigs)} contigs, {genome_len} bp, "
        f"{len(genes)} genes, genomeSAindexNbases={nbases}"
    )

    # Disjoint first-field barcodes between the sublibraries
    half = len(bc_sets["bc1"]) // 2
    seen_r2: set[str] = set()
    seen_cells: set[tuple[str, str, str]] = set()
    for i, (sublib, f1_barcodes) in enumerate(
        [
            ("s1", bc_sets["bc1"][:half]),
            ("s2", bc_sets["bc1"][half:]),
        ]
    ):
        reads = make_reads(
            rng,
            sublib,
            f1_barcodes,
            bc_sets,
            genes,
            args.cells,
            args.reads_per_cell,
            seen_r2,
            seen_cells,
        )
        write_fastqs(reads, outputs[3 + i * 2], outputs[4 + i * 2])
        print(f"# {sublib}: {len(reads)} read pairs, {args.cells} cells")

    # All generated files must stay tiny
    too_big = []
    for path in outputs:
        size = path.stat().st_size
        print(f"# Wrote {path} ({size / 1024:.1f} KB)")
        if size >= MAX_FILE_SIZE:
            too_big.append(path)
    if too_big:
        sys.exit(
            "FATAL: generated file(s) reach 1 MB, the test data must stay "
            "tiny: " + ", ".join(str(path) for path in too_big)
        )


if __name__ == "__main__":
    main()
