from pathlib import Path
from biopipen.utils.misc import run_command, logger

infile = Path({{in.afile | quote}})  # pyright: ignore # noqa: #999
outfile = Path({{in.bfile | quote}})  # pyright: ignore
bedtools: str = {{envs.bedtools | quote}}  # pyright: ignore
window = {{envs.window | repr}}  # pyright: ignore
step = {{envs.step | repr}}  # pyright: ignore
nwin = {{envs.nwin | repr}}  # pyright: ignore
reverse = {{envs.reverse | repr}}  # pyright: ignore
name = {{envs.name | repr}}  # pyright: ignore

if nwin is None and window is None:
    raise ValueError("Either `nwin` or `window` should be provided.")

if nwin is not None and window is not None:
    raise ValueError("Either `nwin` or `window` should be provided, not both.")

# detect if infile is a genome size file or a bed file
with infile.open() as f:
    line = f.readline().strip()
    if len(line.split("\t")) > 2:
        is_bed = True
    else:
        is_bed = False

if is_bed:
    logger.info("BED file is detected as input.")
    cmd = [bedtools, "makewindows", "-b", infile]
else:
    logger.info("Genome size file is detected as input.")
    cmd = [bedtools, "makewindows", "-g", infile]

if nwin:
    cmd.extend(["-n", nwin])
elif step is not None:
    cmd.extend(["-w", window, "-s", step])
else:
    cmd.extend(["-w", window])

if reverse:
    cmd.append("-reverse")

if name != "none":
    cmd.extend(["-name", name])

run_command(cmd, stdout=outfile)
