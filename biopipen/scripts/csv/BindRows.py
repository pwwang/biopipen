import pandas as pd
from biopipen.utils.misc import exec_code

infiles = {{in.infiles | repr}}  # pyright: ignore
outfile = {{in.outfile | repr}}  # pyright: ignore
sep = {{envs.sep | repr}}  # pyright: ignore
header = {{envs.header | repr}}  # pyright: ignore
helper_code = {{envs.helper_code | repr}}  # pyright: ignore
add_filename = {{envs.add_filename | repr}}  # pyright: ignore

exec_code(helper_code, globals(), locals())

if add_filename:
    if not isinstance(add_filename, dict):
        add_filename = {"Filename": add_filename}

    for key, val in add_filename.items():
        add_filename[key] = eval(val)

dfs = []
for infile in infiles:
    df = pd.read_csv(infile, sep=sep, header=header)
    if add_filename:
        for key, val in add_filename.items():
            df[key] = val(infile)
    dfs.append(df)

df = pd.concat(dfs, axis=0)
df.to_csv(outfile, sep=sep, header=header)
