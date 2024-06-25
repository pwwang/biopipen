from biopipen.utils import run_command, dict_to_cli_args

inbed = {{in.inbed | repr}}  # pyright: ignore # noqa: #999
outbed = {{out.outbed | repr}}  # pyright: ignore
envs = {{envs | repr}}  # pyright: ignore
bedtools = envs.pop("bedtools", "bedtools")

envs[""] = [bedtools, "merge"]
envs["i"] = inbed

run_command(dict_to_cli_args(envs, prefix="-"), stdout=outbed)
