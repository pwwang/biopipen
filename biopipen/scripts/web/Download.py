from pathlib import Path
import cmdy

url = {{in.url | repr}}  # pyright: ignore
outfile = Path({{out.outfile | repr}})  # pyright: ignore
tool = {{envs.tool | repr}}  # pyright: ignore
wget = {{envs.wget | repr}}  # pyright: ignore
aria2c = {{envs.aria2c | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore

if tool == "wget":
    args["_"] = url
    args["O"] = outfile
    args["no-directories"] = True
    args["_exe"] = wget
    cmdy.wget(**args).fg()

elif tool == "aria2c":
    args["d"] = outfile.parent
    args["o"] = outfile.name
    args["s"] = ncores
    args["x"] = ncores
    args["_exe"] = aria2c
    args["_"] = url
    cmdy.aria2c(**args).fg()

else: # use python
    import urllib
    urllib.urlretrieve(url, outfile)
