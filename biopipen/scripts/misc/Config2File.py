import json
import rtoml

configstr = {{in.config | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
infmt = {{envs.infmt | quote}}  # pyright: ignore
outfmt = {{envs.outfmt | quote}}  # pyright: ignore

data = rtoml.loads(configstr) if infmt == "toml" else json.loads(configstr)

with open(outfile, "w") as fout:
    if outfmt == "toml":
        rtoml.dump(data, fout)
    else:
        json.dump(data, fout)
