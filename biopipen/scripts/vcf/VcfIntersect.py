from biopipen.utils.misc import run_command, dict_to_cli_args

infile1: str = {{in.infile1 | quote}}  # pyright: ignore  # noqa
infile2 = {{in.infile2 | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
bcftools = {{envs.bcftools | repr}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore
collapse = {{envs.collapse | repr}}  # pyright: ignore

if index:
    gz = True

args = {
    "c": collapse,
    "O": "z" if gz else "v",
    "o": outfile,
    "write": 1,
    "_": ["-n=2", infile1, infile2],
    "": [bcftools, "isec"],
}

run_command(dict_to_cli_args(args, dashify=True), fg=True)

if index:
    run_command([bcftools, "index", "-t", outfile], fg=True)
