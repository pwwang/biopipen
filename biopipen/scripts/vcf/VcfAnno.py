from os import path

from biopipen.utils.misc import run_command, dict_to_cli_args

infile: str = {{in.infile | quote}}  # pyright: ignore  # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore
joboutdir: str = {{job.outdir | quote}}  # pyright: ignore
vcfanno = {{envs.vcfanno | quote}}  # pyright: ignore
ncores: int = {{envs.ncores | repr}}  # pyright: ignore
args: dict = {{envs.args | dict}}  # pyright: ignore

{% set conf = envs.conffile or in.conffile %}  # pyright: ignore  # noqa
{% if conf | isinstance: dict %}  # pyright: ignore  # noqa
conffile = path.join(joboutdir, "config.toml")
conf: str = {{ conf | toml | quote }}  # pyright: ignore  # noqa
with open(conffile, "w") as f:
    f.write(conf)
{% else %}  # pyright: ignore  # noqa
conffile = {{conf | quote}}  # pyright: ignore  # noqa
{% endif %}  # pyright: ignore  # noqa

args["p"] = ncores
args["_"] = [conffile, infile]
args[""] = vcfanno

run_command(dict_to_cli_args(args, dashify=True, prefix="-"), stdout=outfile)
