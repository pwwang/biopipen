import json
import logging
import sys
from pathlib import Path
from prodigy_prot.predict_IC import (  # type: ignore
    Prodigy,
    check_path,
    parse_structure,
)

infile: str = {{in.infile | quote}}  # pyright: ignore # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
distance_cutoff = {{envs.distance_cutoff | float}}  # pyright: ignore
acc_threshold = {{envs.acc_threshold | float}}  # pyright: ignore
temperature = {{envs.temperature | float}}  # pyright: ignore
contact_list = {{envs.contact_list | repr}}  # pyright: ignore
pymol_selection = {{envs.pymol_selection | repr}}  # pyright: ignore
selection = {{envs.selection | repr}}  # pyright: ignore
outtype = {{envs.outtype | repr}}  # pyright: ignore

raw_outfile = Path(outdir) / "_prodigy_raw.txt"
json_outfile = Path(outdir) / "_prodigy.json"
tsv_outfile = Path(outdir) / "_prodigy.tsv"

# log to the raw_outfile
logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
logger = logging.getLogger("Prodigy")

if isinstance(selection, str):
    selection = [selection]

struct_path = check_path(infile)

# parse structure
structure, n_chains, n_res = parse_structure(struct_path)
logger.info(
    "[+] Parsed structure file {0} ({1} chains, {2} residues)".format(
        structure.id, n_chains, n_res
    )
)
prodigy = Prodigy(structure, selection, temperature)
prodigy.predict(distance_cutoff=distance_cutoff, acc_threshold=acc_threshold)
prodigy.print_prediction(outfile=raw_outfile, quiet=False)

# Print out interaction network
if contact_list:
    prodigy.print_contacts(f"{outdir}/prodigy.ic")

# Print out interaction network
if pymol_selection:
    prodigy.print_pymol_script(f"{outdir}/prodigy.pml")

# [+] Reading structure file: <path/to/structure.cif>
# [+] Parsed structure file <structure> (4 chains, 411 residues)
# [+] No. of intermolecular contacts: 191
# [+] No. of charged-charged contacts: 17
# [+] No. of charged-polar contacts: 18
# [+] No. of charged-apolar contacts: 60
# [+] No. of polar-polar contacts: 5
# [+] No. of apolar-polar contacts: 41
# [+] No. of apolar-apolar contacts: 50
# [+] Percentage of apolar NIS residues: 33.90
# [+] Percentage of charged NIS residues: 30.48
# [++] Predicted binding affinity (kcal.mol-1):    -21.3
# [++] Predicted dissociation constant (M) at 25.0˚C:  2.3e-16

output = {}
with open(raw_outfile, "r") as f:
    for line in f:
        if line.startswith("[+"):
            line = line.lstrip("[").lstrip("+").lstrip("]").lstrip()
            if line.startswith("Reading structure file"):
                continue
            if line.startswith("Parsed structure file"):
                continue

            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            if key == "No. of intermolecular contacts":
                output["nIC"] = int(value)
            elif key == "No. of charged-charged contacts":
                output["nCCC"] = int(value)
            elif key == "No. of charged-polar contacts":
                output["nCPC"] = int(value)
            elif key == "No. of charged-apolar contacts":
                output["nCAPC"] = int(value)
            elif key == "No. of polar-polar contacts":
                output["nPPC"] = int(value)
            elif key == "No. of apolar-polar contacts":
                output["nAPPC"] = int(value)
            elif key == "No. of apolar-apolar contacts":
                output["nAPAPC"] = int(value)
            elif key.startswith("Percentage of apolar NIS residues"):
                output["pANISR"] = float(value)
            elif key.startswith("Percentage of charged NIS residues"):
                output["pCNISR"] = float(value)
            elif key.startswith("Predicted binding affinity"):
                output["BindingAffinity"] = float(value)
            elif key.startswith("Predicted dissociation constant"):
                output["DissociationConstant"] = float(value)

with open(json_outfile, "w") as f:
    json.dump(output, f, indent=2)

with open(tsv_outfile, "w") as f:
    f.write("\t".join(output.keys()) + "\n")
    f.write("\t".join(map(str, output.values())) + "\n")

if outtype == "json":
    json_outfile.rename(outfile)
    json_outfile.symlink_to(outfile)
elif outtype == "tsv":
    tsv_outfile.rename(outfile)
    tsv_outfile.symlink_to(outfile)
else:
    raw_outfile.rename(outfile)
    raw_outfile.symlink_to(outfile)
