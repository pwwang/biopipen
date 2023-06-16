import time
import subprocess
from shutil import which
from pathlib import Path, PosixPath  # for as_path

from biopipen.utils.misc import run_command

bamfiles = {{in.bamfiles | repr}}  # pyright: ignore
atfile = {{in.atfile | repr}}  # pyright: ignore

targetfile = {{out.targetfile | repr}}  # pyright: ignore
covfile = {{out.targetfile | as_path | attr: "with_suffix" | call: ".cnn" | repr}}  # pyright: ignore

cnvkit = {{envs.cnvkit | repr}}  # pyright: ignore
samtools = {{envs.samtools | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
ref = {{envs.ref | repr}}  # pyright: ignore
guided = {{envs.guided | repr}}  # pyright: ignore
min_depth = {{envs.min_depth | repr}}  # pyright: ignore
min_gap = {{envs.min_gap | repr}}  # pyright: ignore
min_length = {{envs.min_length | repr}}  # pyright: ignore

params = {}

if guided is None:
    raise ValueError("`envs.guided` is not set.")
elif guided:
    params["t"] = atfile
    params["min-depth"] = min_depth
else:
    params["a"] = atfile
    params["min-gap"] = min_gap
    params["min-length"] = min_length

biopipen_dir = {{biopipen_dir | repr}}  # pyright: ignore

# get the python path from cnvkit.py
cnvkit_path = Path(which(cnvkit))
# Modify cnvkit.py to a unique tmp path, named with timestamp
# to find the python path
tmp_cnvkit_path = Path("/tmp/cnvkit-{}.py".format(time.time()))
with tmp_cnvkit_path.open("w") as f:
    for line in cnvkit_path.read_text().splitlines():
        if line.startswith("if __name__ == "):
            f.write(line + "\n")
            f.write("    import sys\n")
            f.write("    print(sys.executable)\n")
            break
        else:
            f.write(line + "\n")
# make tmp_cnvkit_path executable
tmp_cnvkit_path.chmod(0o755)
# run tmp_cnvkit_path to get the python path
python = subprocess.check_output([tmp_cnvkit_path]).decode("utf-8").strip()

guess_baits = Path(biopipen_dir).joinpath("scripts", "cnvkit", "guess_baits.py")

params.update({
    "o": targetfile,
    "c": covfile,
    "p": ncores,
    "f": ref,
    "s": samtools,
})

cmd = [python, guess_baits]
for k, v in params.items():
    if len(k) == 1:
        cmd.append("-{}".format(k))
        cmd.append(v)
    else:
        cmd.append("--{}".format(k))
        cmd.append(v)

cmd.extend(bamfiles)

run_command(cmd, fg=True)
