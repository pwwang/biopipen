from pathlib import Path
from biopipen.utils.misc import logger, run_command, dict_to_cli_args

nsnps = {{in.nsnps | repr}}  # pyright: ignore
ncases = {{in.ncases | repr}}  # pyright: ignore
nctrls = {{in.nctrls | repr}}  # pyright: ignore
outdir = {{out.outdir | repr}}  # pyright: ignore
plink = {{envs.plink | repr}}  # pyright: ignore
seed = {{envs.seed | repr}}  # pyright: ignore
label = {{envs.label | repr}}  # pyright: ignore
prevalence = {{envs.prevalence | repr}}  # pyright: ignore
minfreq = {{envs.minfreq | repr}}  # pyright: ignore
maxfreq = {{envs.maxfreq | repr}}  # pyright: ignore
hetodds = {{envs.hetodds | repr}}  # pyright: ignore
homodds = {{envs.homodds | repr}}  # pyright: ignore
missing = {{envs.missing | repr}}  # pyright: ignore
args = {{envs.args | repr}}  # pyright: ignore

logger.info("Generating parameters file")
params_file = Path(outdir) / "params.txt"
params_file.write_text(
    f"{nsnps}\t{label}\t{minfreq}\t{maxfreq}\t{hetodds}\t{homodds}\n"
)

if seed is not None:
    args["seed"] = seed

args["simulate"] = params_file
args["out"] = Path(outdir) / "sim_snps"
args["simulate-ncases"] = ncases
args["simulate-ncontrols"] = nctrls
args["simulate-prevalence"] = prevalence
args["simulate-missing"] = missing

cmd = [plink] + dict_to_cli_args(args)

logger.info("Running PLINK simulation ...")
run_command(cmd, fg=True)
