from biopipen.utils.misc import run_command, dict_to_cli_args

infile = {{in.infile | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
bcftools = {{envs.bcftools | quote}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore
tmpdir = {{envs.tmpdir | quote}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore

args[""] = bcftools
args["_"] = infile
args["o"] = outfile
args["O"] = "z" if gz or index else "v"

run_command(dict_to_cli_args(args, dashify=True), fg=True)

if index:
    run_command([bcftools, "index", outfile], fg=True)
