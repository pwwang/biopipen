from biopipen.utils.misc import run_command, dict_to_cli_args

inbed = {{in.inbed | quote}}  # pyright: ignore # noqa: #999
outbed = {{out.outbed | quote}}  # pyright: ignore
envs: dict = {{envs | dict}}  # pyright: ignore
bedtools = envs.pop("bedtools", "bedtools")

envs[""] = [bedtools, "merge"]
envs["i"] = inbed

run_command(dict_to_cli_args(envs, prefix="-"), stdout=outbed)
