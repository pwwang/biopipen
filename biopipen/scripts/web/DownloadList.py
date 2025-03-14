from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args

urlfile: str = {{in.urlfile | quote}}  # pyright: ignore  # noqa
outdir = Path({{out.outdir | quote}})  # pyright: ignore
tool = {{envs.tool | repr}}  # pyright: ignore
wget = {{envs.wget | repr}}  # pyright: ignore
aria2c = {{envs.aria2c | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
args: dict = {{envs.args | repr}}  # pyright: ignore

if tool == "wget":
    args["i"] = urlfile
    args["P"] = outdir
    args["no-directories"] = True
    args[""] = wget
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

elif tool == "aria2c":
    args["d"] = outdir
    args["s"] = ncores
    args["x"] = ncores
    args[""] = aria2c
    args["i"] = urlfile
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

else: # use python
    from urllib.request import urlretrieve
    from urllib.parse import urlparse
    with open(urlfile, "r") as furl:
        for i, url in enumerate(furl.readlines()):
            parsed = urlparse(url)
            path = Path(parsed.path)
            urlretrieve(url, f"{path.stem}-{i}{path.suffix}")
