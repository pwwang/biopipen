from pathlib import Path
from biopipen.utils.misc import logger, run_command, dict_to_cli_args

nsnps = {{in.nsnps | repr}}  # pyright: ignore
ncases = {{in.ncases | repr}}  # pyright: ignore
nctrls = {{in.nctrls | repr}}  # pyright: ignore
outdir = {{out.outdir | repr}}  # pyright: ignore
gtmatfile = {{out.gtmat | repr}}  # pyright: ignore
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
transpose_gtmat = {{envs.transpose_gtmat | repr}}  # pyright: ignore
sample_prefix = {{envs.sample_prefix | repr}}  # pyright: ignore

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

# Transpose the genotype matrix
# CHR	SNP	(C)M	POS	COUNTED	ALT	per0_per0	per1_per1	per2_per2
# 1	SNP_0	0	1	D	d	1	0	1
# 1	SNP_1	0	2	d	D	0	1	0
# 1	SNP_2	0	3	d	D	0	0	0
# 1	SNP_3	0	4	d	D	0	0	0
# 1	SNP_4	0	5	D	d	1	2	1
cmd = [
    plink,
    "--recode",
    "A" if transpose_gtmat else "A-transpose",
    "tab",
    "--bfile",
    args["out"],
    "--out",
    gtmatfile + ".plink.recoded",
]
logger.info("Recoding into genotype matrix ...")
run_command(cmd, fg=True)

logger.info("Saving genotype matrix ...")
## transpose_gtmat = False
# SNP_COUNTED	per0_per0	per1_per1	per2_per2
# SNP_0_D	1	0	1
# SNP_1_d	0	1	0
# SNP_2_d	0	0	0
# SNP_3_d	0	0	0
# SNP_4_D	1	2	1
## transpose_gtmat = True
# FID_IID SNP_0_D SNP_1_D SNP_2_D
# per0_per0 0 1 1
# per1_per1 0 2 0
# per2_per2 0 0 0
# per3_per3 1 1 0
# per4_per4 0 0 0
if transpose_gtmat:
    cmd = f"cut -f1,2,7- {gtmatfile}.plink.recoded.raw | sed 's/\\t/_/'"
else:
    cmd = f"cut -f2,5,7- {gtmatfile}.plink.recoded.traw | sed 's/\\t/_/'"

if sample_prefix:
    cmd = f"{cmd} | sed 's/per[0-9]\\+_per/{sample_prefix}/g'"

cmd = f"{cmd} > {gtmatfile}"

run_command(cmd, fg=True)
