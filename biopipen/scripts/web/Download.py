from pathlib import Path
import cmdy

url = {{in.url | repr}}
outfile = Path({{out.outfile | repr}})
tool = {{envs.tool | repr}}
wget = {{envs.wget | repr}}
aria2c = {{envs.aria2c | repr}}
ncores = {{envs.ncores | repr}}
args = {{envs.args | repr}}

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
