from pathlib import Path
from shutil import which
from diot import Diot  # noqa: F401
from biopipen.utils.misc import run_command, dict_to_cli_args

infile: str = {{in.infile | quote}}  # pyright: ignore # noqa
outfile: str = {{out.outfile | quote}}  # pyright: ignore
envs: dict = {{envs | repr}}  # pyright: ignore
tool: str = envs.pop("tool", "maxit")
maxit: str = envs.pop("maxit", "maxit")
beem = envs.pop("beem", "BeEM")

if tool == "maxit":
    maxit_found = which(maxit)
    if not maxit_found:
        raise ValueError(f"maxit executable not found: {maxit}")

    maxit_exe = Path(maxit_found).expanduser().resolve()
    rcsbroot = maxit_exe.parent.parent
    envs["input"] = infile
    envs["output"] = outfile
    envs["o"] = 2
    envs["log"] = Path(outfile).with_suffix(".log")
    run_command([maxit, *dict_to_cli_args(envs, prefix="-")], fg=True, env={"RCSBROOT": rcsbroot})

else:
    outfile: Path = Path(outfile)  # type: ignore
    envs["_"] = infile
    envs["p"] = outfile.parent.joinpath(outfile.stem)
    envs["outfmt"] = 3
    args = dict_to_cli_args(envs, prefix="-", sep="=")

    run_command([beem, *args], fg=True)
