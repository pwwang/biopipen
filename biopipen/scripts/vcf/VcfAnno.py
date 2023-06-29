from os import path

from biopipen.utils.misc import run_command, dict_to_cli_args

infile = {{in.infile | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
joboutdir = {{job.outdir | quote}}  # pyright: ignore
vcfanno = {{envs.vcfanno | quote}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore

{% set conf = envs.conffile or in.conffile %}
{% if conf | isinstance: dict %}
conffile = path.join(joboutdir, "config.toml")
conf = {{ conf | toml | quote }}
with open(conffile, "w") as f:
    f.write(conf)
{% else %}
conffile = {{conf | quote}}
{% endif %}

args["p"] = ncores
args["_"] = [conffile, infile]
args[""] = vcfanno

run_command(dict_to_cli_args(args, dashify=True, prefix="-"), stdout=outfile)
