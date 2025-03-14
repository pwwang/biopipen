import json
import rtoml

configstr: str = {{in.config | quote}}  # pyright: ignore  # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore
infmt = {{envs.infmt | quote}}  # pyright: ignore
outfmt = {{envs.outfmt | quote}}  # pyright: ignore

data = rtoml.loads(configstr) if infmt == "toml" else json.loads(configstr)

with open(outfile, "w") as fout:
    if outfmt == "toml":
        rtoml.dump(data, fout)
    else:
        json.dump(data, fout)
