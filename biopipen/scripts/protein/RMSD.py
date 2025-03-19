from pathlib import Path
from shutil import which
from diot import Diot  # noqa: F401
from biopipen.utils.misc import run_command, dict_to_cli_args

infile1: str = {{in.infile1 | quote}}  # pyright: ignore # noqa
infile2: str = {{in.infile2 | quote}}  # pyright: ignore # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore # noqa
outdir: str = {{job.outdir | quote}}  # pyright: ignore # noqa
envs: dict = {{envs | repr}}  # pyright: ignore # noqa
conv_tool = envs.pop("conv_tool", "maxit")
maxit = envs.pop("maxit", "maxit")
beem = envs.pop("beem", "BeEM")
ca_only = envs.pop("ca_only", False)
# aa20_only = envs.pop("aa20_only", False)
duel = envs.pop("duel", "keep")
calculate_rmsd = envs.pop("calculate_rmsd", "calculate_rmsd")


def cif_to_pdb(cif_file, pdb_file:Path):
    if conv_tool == "maxit":
        maxit_bin = Path(which(maxit)).resolve()
        rcsbroot = Path(maxit_bin).parent.parent
        args = {"input": cif_file, "output": pdb_file, "o": 2, "log": pdb_file.with_suffix(".log")}
        run_command([maxit, *dict_to_cli_args(args, prefix="-")], fg=True, env={"RCSBROOT": rcsbroot})
    else:
        args = {"_": cif_file, "p": pdb_file.parent.joinpath(pdb_file.stem)}
        args = dict_to_cli_args(args, prefix="-", sep="=")
        run_command([beem, *args], fg=True)


def pdb_to_ca_pdb(pdb_file: Path, ca_pdb_file: Path):
    """Extract C-alpha atoms from a PDB file and still keep the original order and metadata."""
    with open(pdb_file, "r") as f, open(ca_pdb_file, "w") as fw:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                fw.write(line)


# def pdb_to_aa20_pdb(pdb_file: Path, aa20_pdb_file: Path):
#     """Extract the 20 amino acids from a PDB file and still keep the original order and metadata."""
#     with open(pdb_file, "r") as f, open(aa20_pdb_file, "w") as fw:
#         for line in f:
#             if line.startswith("ATOM") and line[17:20].strip() in (
#                 "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
#                 "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
#                 "THR", "TRP", "TYR", "VAL",
#             ):
#                 fw.write(line)


def deduel_pdb(pdb_file: Path, deduel_pdb_file: Path):
    """Remove/Handle the duel atoms in a PDB file."""
    def is_duel(atom1, atom2):
        #           1         2
        # 01234567890123456789012345
        # ATOM    913  CA ATYR A 113
        # ATOM    914  CA BTYR A 113
        # The key should be "ATOM|CA |TYR| A| 113"
        return (
            atom1[:4] == atom2[:4] and
            atom1[12:16] == atom2[12:16] and
            atom1[17:20] == atom2[17:20] and
            atom1[21] == atom2[21] and
            atom1[22:26] == atom2[22:26] and
            atom1[16] != atom2[16]
        )

    def clean_atom(atom):
        return atom[:16] + " " + atom[17:]

    last_atom = ""
    with open(pdb_file, "r") as f, open(deduel_pdb_file, "w") as fw:
        for line in f:
            if not line.startswith("ATOM"):
                fw.write(line)
                continue
            if not is_duel(last_atom, line):
                if last_atom:
                    fw.write(clean_atom(last_atom))
                last_atom = line
            # is duel
            elif duel == "keep":
                fw.write(clean_atom(last_atom))
                fw.write(clean_atom(line))
                last_atom = ""
            elif duel == "keep_first":
                fw.write(clean_atom(last_atom))
                last_atom = ""
            elif duel == "keep_last":
                fw.write(clean_atom(line))
                last_atom = ""
            elif duel == "average":
                # Average the coordinates
                x1 = float(last_atom[30:38])
                y1 = float(last_atom[38:46])
                z1 = float(last_atom[46:54])
                x2 = float(line[30:38])
                y2 = float(line[38:46])
                z2 = float(line[46:54])
                x = (x1 + x2) / 2.0
                y = (y1 + y2) / 2.0
                z = (z1 + z2) / 2.0
                fw.write(clean_atom(last_atom[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + last_atom[54:]))
                last_atom = ""

        if last_atom:
            fw.write(last_atom)


def index_of(lst, item) -> int:
    try:
        return lst.index(item)
    except ValueError:
        return -1


if infile1.endswith(".cif"):
    pdb1 = Path(outdir) / f"{Path(infile1).stem}.pdb"
    cif_to_pdb(infile1, pdb1)
    infile1 = pdb1  # type: ignore

if infile2.endswith(".cif"):
    pdb2 = Path(outdir) / f"{Path(infile2).stem}.pdb"
    cif_to_pdb(infile2, pdb2)
    infile2 = pdb2  # type: ignore

if ca_only:
    ca_pdb1 = Path(outdir) / f"{Path(infile1).stem}.ca.pdb"
    pdb_to_ca_pdb(infile1, ca_pdb1) # type: ignore
    infile1 = ca_pdb1  # type: ignore

    ca_pdb2 = Path(outdir) / f"{Path(infile2).stem}.ca.pdb"
    pdb_to_ca_pdb(infile2, ca_pdb2) # type: ignore
    infile2 = ca_pdb2  # type: ignore

# if aa20_only:
#     aa20_pdb1 = Path(outdir) / f"{Path(infile1).stem}.aa20.pdb"
#     pdb_to_aa20_pdb(infile1, aa20_pdb1) # type: ignore
#     infile1 = aa20_pdb1  # type: ignore

#     aa20_pdb2 = Path(outdir) / f"{Path(infile2).stem}.aa20.pdb"
#     pdb_to_aa20_pdb(infile2, aa20_pdb2) # type: ignore
#     infile2 = aa20_pdb2  # type: ignore

if duel != "keep":
    deduel_pdb1 = Path(outdir) / f"{Path(infile1).stem}.deduel.pdb"
    deduel_pdb(infile1, deduel_pdb1) # type: ignore
    infile1 = deduel_pdb1  # type: ignore

    deduel_pdb2 = Path(outdir) / f"{Path(infile2).stem}.deduel.pdb"
    deduel_pdb(infile2, deduel_pdb2) # type: ignore
    infile2 = deduel_pdb2  # type: ignore

envs["_"] = [infile1, infile2]
envs = dict_to_cli_args(envs, dashify=True)

idx_ur = index_of(envs, "--ur")
if idx_ur != -1:
    envs[idx_ur] = "-ur"

idx_urks = index_of(envs, "--urks")
if idx_urks != -1:
    envs[idx_urks] = "-urks"

idx_nh = index_of(envs, "--nh")
if idx_nh != -1:
    envs[idx_nh] = "-nh"

out: str = run_command([calculate_rmsd, *envs], stdout="return")  # type: ignore
out = out.strip()

try:
    float(out)
except (ValueError, TypeError):
    raise ValueError(out)

Path(outfile).write_text(out)
