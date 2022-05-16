import json
import rtoml

configstr = {{in.config | repr}}
outfile = {{out.outfile | quote}}
infmt = {{envs.infmt | quote}}
outfmt = {{envs.outfmt | quote}}

data = rtoml.loads(configstr) if infmt == "toml" else json.loads(configstr)

with open(outfile, "w") as fout:
    if outfmt == "toml":
        rtoml.dump(data, fout)
    else:
        json.dump(data, fout)
