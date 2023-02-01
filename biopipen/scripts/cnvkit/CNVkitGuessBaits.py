from pathlib import Path, PosixPath  # for as_path

import cmdy

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
python = Path(
    cmdy.which(cnvkit).stdout.strip()
).resolve().read_text().splitlines()[0][2:].strip()
guess_baits = Path(biopipen_dir).joinpath("scripts", "cnvkit", "guess_baits.py")

params.update({
    "o": targetfile,
    "c": covfile,
    "p": ncores,
    "f": ref,
    "s": samtools,
    "_": bamfiles,
    "_exe": python
})

cmd = cmdy.python(guess_baits, **params).hold()
print("Running command:")
print(cmd.strcmd)
cmd.fg().run()
