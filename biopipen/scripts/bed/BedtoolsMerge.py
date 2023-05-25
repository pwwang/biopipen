import shlex
from biopipen.utils import run_command

inbed = {{in.inbed | repr}}  # pyright: ignore
outbed = {{out.outbed | repr}}  # pyright: ignore
envs = {{envs | repr}}  # pyright: ignore
bedtools = envs.pop("bedtools", "bedtools")

cmd = f"{shlex.quote(bedtools)} merge -i {shlex.quote(inbed)}"
for k, v in envs.items():
    if v is True:
        cmd += f" -{k}"
    elif v is not False:
        cmd += f" -{k} {shlex.quote(v)}"

cmd += f" > {shlex.quote(outbed)}"

run_command(cmd, fg=True, print_command=True)
