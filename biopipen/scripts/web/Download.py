from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args

url = {{in.url | quote}}  # pyright: ignore # noqa
outfile = Path({{out.outfile | quote}})  # pyright: ignore
tool = {{envs.tool | repr}}  # pyright: ignore
wget = {{envs.wget | repr}}  # pyright: ignore
aria2c = {{envs.aria2c | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
args: dict = {{envs.args | dict}}  # pyright: ignore

if tool == "wget":
    args["_"] = url
    args["O"] = outfile
    args["no-directories"] = True
    args[""] = wget
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

elif tool == "aria2c":
    args["d"] = outfile.parent
    args["o"] = outfile.name
    args["s"] = ncores
    args["x"] = ncores
    args[""] = aria2c
    args["_"] = url
    run_command(dict_to_cli_args(args, dashify=True), fg=True)

else: # use python
    import urllib

    try:
        urllib.urlretrieve(url, outfile)  # type: ignore
    except AttributeError:
        urllib.request.urlretrieve(url, outfile)  # type: ignore
