# """
# LICENSE

# GNU General Public License v2.0

# The code is based on the script from:
# https://github.com/kad-ecoli/pdb2fasta/blob/master/pdb2fasta.py

# The original code is licensed under GNU General Public License v2.0.
# The original code is modified by biopipen developers to fit the biopipen.
# """
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path

infile: str = {{in.infile | quote}}  # pyright: ignore # noqa: E999
outfile: str = {{out.outfile | quote}}  # pyright: ignore
chains: str | list | None = {{envs.chains | repr}}  # pyright: ignore
wrap: int = {{envs.wrap | repr}}  # pyright: ignore

if isinstance(chains, str):
    chains = [chain.strip() for chain in chains.split(",")]

aa3to1 = {
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
   'MSE':'M',
}

ca_pattern = re.compile(
    r"^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])"  # noqa: W605
)

filename = Path(infile).stem
chain_dict = defaultdict(str)

with open(infile, 'r') as fp:
    for line in fp:
        if line.startswith("ENDMDL"):
            break

        match_list = ca_pattern.findall(line)
        if match_list:
            resn = match_list[0][0] + match_list[0][2]
            chain = match_list[0][1] + match_list[0][3]
            if chains is None or chain in chains:
                chain_dict[chain] += aa3to1[resn]

with open(outfile, 'w') as fp:
    for chain in chain_dict:
        fp.write(f">{filename}:{chain}\n")
        sequence = chain_dict[chain]
        if wrap > 0:
            for i in range(0, len(sequence), 80):
                fp.write(sequence[i:i+80] + "\n")
        else:
            fp.write(sequence + "\n")
