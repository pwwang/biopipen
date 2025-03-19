"""Script for regulatory.MotifScan"""
import re

# Paths may be passed in args or to motifdb
from pathlib import PosixPath  # noqa: F401
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

motiffile: str = {{in.motiffile | quote}}  # pyright: ignore # noqa: #999
seqfile: str = {{in.seqfile | quote}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore

tool = {{envs.tool | repr}}  # pyright: ignore
fimo = {{envs.fimo | repr}}  # pyright: ignore
motif_col: str | int = {{envs.motif_col | repr}}  # pyright: ignore
regulator_col: str | int = {{envs.regulator_col | repr}}  # pyright: ignore
notfound = {{envs.notfound | repr}}  # pyright: ignore
motifdb: str | None = {{envs.motifdb | repr}}  # pyright: ignore
cutoff = {{envs.cutoff | repr}}  # pyright: ignore
q = {{envs.q | repr}}  # pyright: ignore
q_cutoff = {{envs.q_cutoff | repr}}  # pyright: ignore
args: dict = {{envs.args | dict | repr}}  # pyright: ignore

# Check if the tool is supported
if tool != "fimo":
    raise ValueError(f"Unsupported tool: {tool}, currently only fimo is supported")

# Check if the motif database is provided
if motifdb is None:
    raise ValueError("The motif database is required")

# Check if the motif file exists
if not motiffile:
    raise FileNotFoundError(f"Motif file in.motiffile must be provided")

# Check if the sequence file exists
if not seqfile:
    raise FileNotFoundError(f"Sequence file in.seqfile must be provided")

# Normalize motif_col and regulator_col into 0-based indexes
if isinstance(motif_col, str) or isinstance(regulator_col, str):
    with open(motiffile, "r") as f:
        header = f.readline().strip().split("\t")
        if isinstance(motif_col, str):
            motif_col: int = header.index(motif_col) + 1
        if isinstance(regulator_col, str):
            regulator_col = header.index(regulator_col) + 1
if isinstance(motif_col, int):
    motif_col -= 1
if isinstance(regulator_col, int):
    regulator_col -= 1

# Check if motif names exist in the database
with open(motiffile, "r") as f:
    motif_names = set(
        line.strip().split("\t")[motif_col]
        for i, line in enumerate(f)
        if i > 0  # skip header
    )

with open(motifdb, "r") as f:
    motif_db_names = set(
        line[6:].strip()
        for line in f
        if line.startswith("MOTIF")
    )

if notfound == "error":
    notfound_motifs = motif_names - motif_db_names
    if notfound_motifs:
        raise ValueError(f"Motifs not found in the database: {notfound_motifs}")

# Make a new motif database with only the motifs in the motiffile
motif_names = motif_names & motif_db_names
motifdb_filtered = f"{outdir}/motif_db.txt"
with open(motifdb, "r") as f, open(motifdb_filtered, "w") as f_out:
    should_write = True
    for line in f:
        if line.startswith("MOTIF"):
            motif_name = line[6:].strip()
            if motif_name in motif_names:
                should_write = True
            else:
                should_write = False

        if should_write:
            f_out.write(line)
        else:
            continue

# Now run fimo
args[""] = fimo
args["oc"] = f"{outdir}"
args["thresh"] = cutoff
args["qv_thresh"] = q_cutoff
args["no_qvalue"] = not q
args["no-pgc"] = True
args["_"] = [motifdb_filtered, seqfile]

logger.info("Running fimo ...")
run_command(dict_to_cli_args(args, dashify=True), fg=True)

logger.info("Adding additional information to the output ...")
# Get the motif to regulator mapping
motif_regulator_map = {}
if regulator_col is not None:
    with open(motiffile, "r") as f:
        next(f)  # skip header
        for line in f:
            line = line.strip().split("\t")
            motif_name = line[motif_col]
            regulator = line[regulator_col]
            motif_regulator_map[motif_name] = regulator

# Get the sequence name information
seqnames = {}
seqcoords = {}
with open(seqfile, "r") as f:
    for line in f:
        if not line.startswith(">"):
            continue

        seqname = line[1:].strip()
        match = re.match(r"^(.+)::((?:chr)?\d+):(\d+)-(\d+).*$", seqname)
        if not match:
            seqnames[seqname] = seqname
            seqcoords[seqname] = None
        else:
            sname, chrom, start, end = match.groups()
            seqnames[seqname] = sname
            seqcoords[seqname] = (chrom, int(start), int(end))

# Add additional information to the output
with open(f"{outdir}/fimo.tsv", "r") as f, open(f"{outdir}/fimo_output.txt", "w") as f_out:
    header = f.readline().strip().split("\t")
    f_out.write(
        "\t".join(header + ["regulator", "seqname", "seqstart", "seqstop"]) + "\n"
    )
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        line = line.split("\t")
        motif_name = line[0]
        sequence_name = line[2]
        start = int(line[3])
        stop = int(line[4])
        regulator = motif_regulator_map.get(motif_name, motif_name)
        seqname = seqnames.get(sequence_name, "NA")
        seqcoord = seqcoords.get(sequence_name)
        if not seqcoord:
            seqstart = "NA"
            seqstop = "NA"
        else:
            seqstart = start + seqcoord[1] - 1
            seqstop = stop + seqcoord[2] - 1

        f_out.write(
            "\t".join(line + [regulator, seqname, str(seqstart), str(seqstop)]) + "\n"
        )
