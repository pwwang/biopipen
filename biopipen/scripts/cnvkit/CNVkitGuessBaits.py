import time
import subprocess
from shutil import which
from pathlib import Path, PosixPath  # for as_path

from biopipen.utils.misc import run_command, dict_to_cli_args

bamfiles = {{in.bamfiles | each: str | repr}}  # pyright: ignore  # noqa
atfile = {{in.atfile | quote}}  # pyright: ignore

targetfile = {{out.targetfile | quote}}  # pyright: ignore
covfile = {{out.targetfile | as_path | attr: "with_suffix" | call: ".cnn" | repr}}  # pyright: ignore

cnvkit: str = {{envs.cnvkit | repr}}  # pyright: ignore
samtools = {{envs.samtools | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
ref: str = {{envs.ref | repr}}  # pyright: ignore
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

biopipen_dir: str = {{biopipen_dir | quote}}  # pyright: ignore

# get the python path from cnvkit.py
cnvkit_found = which(cnvkit)
if cnvkit_found is None:
    raise ValueError(f"cnvkit executable not found: {cnvkit}")

cnvkit_path = Path(cnvkit_found)
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
    "": [python, guess_baits],
    "o": targetfile,
    "c": covfile,
    "p": ncores,
    "f": Path(ref).expanduser(),
    "s": samtools,
    "_": bamfiles,
})

run_command(dict_to_cli_args(params), fg=True)
