from __future__ import annotations

import os
import sys
from pathlib import Path

import rtoml
import pandas as pd
from tempfile import gettempdir
from biopipen.utils.misc import logger, run_command

configfile: str = {{in.configfile | quote}}  # pyright: ignore  # noqa
outdir = Path({{out.outdir | quote}})  # pyright: ignore
envs: dict = {{envs | dict | repr}}  # pyright: ignore
python: str | list[str] = sys.executable

args = envs.copy()
config = rtoml.load(Path(configfile))
args.update(config)
model_name = args.pop("model_name")
model_file = Path(args.pop("model_file"))
data_dir = args.pop("data_dir", None)
tcrdock: Path | str | None = args.pop("tcrdock", None)
tmpdir: str = args.pop("tmpdir", gettempdir())
python = args.pop("python", python)

if not isinstance(python, (list, tuple)):
    python = [python]

if not data_dir:
    raise ValueError("`envs.data_dir` is required")

if not tcrdock:
    logger.info("- `envs.tcrdock` is not provided, cloning the repository ... ")
    repo_url = "https://github.com/phbradley/TCRdock"
    commit_id = "c5a7af42eeb0c2a4492a4d4fe803f1f9aafb6193"
    branch = "main"

    from git import Repo
    repo = Repo.clone_from(repo_url, tmpdir, branch=branch, no_checkout=True)
    repo.git.checkout(commit_id)
    tcrdock = Path(tmpdir) / "TCRdock"

    logger.info("- Running download_blast.py ...")
    cmd = [
        *python,
        tcrdock / "download_blast.py",
    ]
    run_command(cmd, fg=True, cwd=str(tcrdock))

tcrdock = str(tcrdock)

if not model_file.is_absolute():
    model_file = Path(data_dir) / "params" / model_file

os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'

logger.info("- Composing targets file ... ")
targets_file = outdir / "user_targets.tsv"
targets = pd.DataFrame(
    [
        dict(
            organism=args['organism'],
            mhc_class=args['mhc_class'],
            mhc=args['mhc'],
            peptide=args['peptide'],
            va=args['va'],
            ja=args['ja'],
            cdr3a=args['cdr3a'],
            vb=args['vb'],
            jb=args['jb'],
            cdr3b=args['cdr3b'],
        )
    ]
)
targets.to_csv(targets_file, sep="\t", index=False)

logger.info("- Generating inputs for AlphaFold modeling ... ")
cmd = [
    *python,
    tcrdock + "/setup_for_alphafold.py",
    "--targets_tsvfile", targets_file,
    "--output_dir", outdir / "user_output",
    "--new_docking",
]
run_command(cmd, fg=True)

logger.info("- Running AlphaFold modeling ... ")
cmd = [
    *python,
    tcrdock + "/run_prediction.py",
    "--verbose",
    "--targets", outdir / "user_output/targets.tsv",
    "--outfile_prefix", f"{outdir}/{args['peptide']}",
    "--model_names", model_name,
    "--data_dir", data_dir,
    "--model_params_files", model_file,
]
run_command(cmd, fg=True, env={"XLA_FLAGS": "--xla_gpu_force_compilation_parallelism=1"})

logger.info("- Calculating the PAE ... ")
cmd = [
    *python,
    tcrdock + "/add_pmhc_tcr_pae_to_tsvfile.py",
    "--infile", f"{outdir}/{args['peptide']}_final.tsv",
    "--outfile", f"{outdir}/{args['peptide']}_w_pae.tsv",
]

run_command(cmd, fg=True)
