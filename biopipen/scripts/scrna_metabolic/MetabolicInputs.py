from pathlib import Path

import rtoml
from diot import Diot

in_metafile = Path({{in.metafile | quote}})        # pyright: ignore
in_gmtfile = Path({{in.gmtfile | quote}})          # pyright: ignore
in_config = {{in.config | repr}}                   # pyright: ignore
out_metafile = Path({{out.metafile | quote}})      # pyright: ignore
out_gmtfile = Path({{out.gmtfile | quote}})        # pyright: ignore
out_configfile = Path({{out.configfile | quote}})  # pyright: ignore

# make symbolic link of in.metafile to out.metafile
out_metafile.symlink_to(in_metafile)

# make symbolic link of in.gmtfile to out.gmtfile
out_gmtfile.symlink_to(in_gmtfile)

# Check the config and save as toml file
if isinstance(in_config, str):
    if "\n" not in in_config and Path(in_config).is_file():
        Path(in_config).symlink_to(out_configfile)
    else:  # A toml string
        with out_configfile.open("w") as f:
            f.write(in_config)
else:  # A python dictionary
    with out_configfile.open("w") as f:
        rtoml.dump(in_config, f)
