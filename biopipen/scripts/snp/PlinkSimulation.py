from pathlib import Path
from multiprocessing import Pool
from slugify import slugify
from simpleconf import Config
from biopipen.utils.misc import logger, run_command, dict_to_cli_args

configfile: str = {{in.configfile | quote}}  # pyright: ignore # noqa: E999
outdir: str = {{out.outdir | quote}}  # pyright: ignore
gtmatfile: str = {{out.gtmat | quote}}  # pyright: ignore
config = Config.load(configfile)

default_nsnps = {{envs.nsnps | repr}}  # pyright: ignore
default_ncases = {{envs.ncases | repr}}  # pyright: ignore
default_nctrls = {{envs.nctrls | repr}}  # pyright: ignore
default_plink = {{envs.plink | repr}}  # pyright: ignore
default_seed = {{envs.seed | repr}}  # pyright: ignore
default_label = {{envs.label | repr}}  # pyright: ignore
default_prevalence = {{envs.prevalence | repr}}  # pyright: ignore
default_minfreq = {{envs.minfreq | repr}}  # pyright: ignore
default_maxfreq = {{envs.maxfreq | repr}}  # pyright: ignore
default_hetodds = {{envs.hetodds | repr}}  # pyright: ignore
default_homodds = {{envs.homodds | repr}}  # pyright: ignore
default_missing = {{envs.missing | repr}}  # pyright: ignore
default_args: dict = {{envs.args | repr}}  # pyright: ignore
default_transpose_gtmat = {{envs.transpose_gtmat | repr}}  # pyright: ignore
default_sample_prefix = {{envs.sample_prefix | repr}}  # pyright: ignore

defaults = {
    "nsnps": default_nsnps,
    "ncases": default_ncases,
    "nctrls": default_nctrls,
    "plink": default_plink,
    "seed": default_seed,
    "label": default_label,
    "prevalence": default_prevalence,
    "minfreq": default_minfreq,
    "maxfreq": default_maxfreq,
    "hetodds": default_hetodds,
    "homodds": default_homodds,
    "missing": default_missing,
    # "args": default_args,
    "transpose_gtmat": default_transpose_gtmat,
    "sample_prefix": default_sample_prefix,
}

def do_one_simulation(confitems):
    args = default_args.copy()
    args.update(confitems.pop("args", {}))
    confs = defaults.copy()
    confs.update(confitems)
    transpose_gtmat = confs.pop("transpose_gtmat")
    sample_prefix = confs.pop("sample_prefix")


    logger.debug("  Generating parameters file")
    params_file = Path(outdir) / "params.txt"
    params_file.write_text(
        f"{confs['nsnps']}\t{confs['label']}\t{confs['minfreq']}\t"
        f"{confs['maxfreq']}\t{confs['hetodds']}\t{confs['homodds']}\n"
    )

    if confs.get('seed') is not None:
        args["seed"] = confs['seed']

    args["simulate"] = params_file
    args["out"] = Path(outdir) / "sim_snps"
    args["simulate-ncases"] = confs['ncases']
    args["simulate-ncontrols"] = confs['nctrls']
    args["simulate-prevalence"] = confs['prevalence']
    args["simulate-missing"] = confs['missing']

    cmd = [confs['plink']] + dict_to_cli_args(args)

    logger.debug("  Running PLINK simulation ...")
    run_command(cmd, fg=True)

    # Transpose the genotype matrix
    # CHR	SNP	(C)M	POS	COUNTED	ALT	per0_per0	per1_per1	per2_per2
    # 1	SNP_0	0	1	D	d	1	0	1
    # 1	SNP_1	0	2	d	D	0	1	0
    # 1	SNP_2	0	3	d	D	0	0	0
    # 1	SNP_3	0	4	d	D	0	0	0
    # 1	SNP_4	0	5	D	d	1	2	1
    cmd = [
        confs['plink'],
        "--recode",
        "A" if transpose_gtmat else "A-transpose",
        "tab",
        "--bfile",
        args["out"],
        "--out",
        gtmatfile + ".plink.recoded",
    ]
    logger.debug("- Recoding into genotype matrix ...")
    run_command(cmd, fg=True)

    logger.debug("  Saving genotype matrix ...")
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


do_one_simulation(config)
