from pathlib import Path
import cmdy

urlfile = {{in.urlfile | repr}}  # pyright: ignore
outdir = Path({{out.outdir | repr}})  # pyright: ignore
tool = {{envs.tool | repr}}  # pyright: ignore
wget = {{envs.wget | repr}}  # pyright: ignore
aria2c = {{envs.aria2c | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore

if tool == "wget":
    args["i"] = urlfile
    args["P"] = outdir
    args["no-directories"] = True
    args["_exe"] = wget
    cmdy.wget(**args).fg()

elif tool == "aria2c":
    args["d"] = outdir
    args["s"] = ncores
    args["x"] = ncores
    args["_exe"] = aria2c
    args["i"] = urlfile
    cmdy.aria2c(**args).fg()

else: # use python
    import urllib
    from urllib.parse import urlparse
    with open(urlfile, "r") as furl:
        for i, url in enumerate(furl.readlines()):
            parsed = urlparse(url)
            path = Path(parsed.path)
            urllib.urlretrieve(url, f"{path.stem}-{i}{path.suffix}")
