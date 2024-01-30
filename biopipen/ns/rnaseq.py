"""RNA-seq data analysis"""

from ..core.proc import Proc
from ..core.config import config


class UnitConversion(Proc):
    """Convert expression value units back and forth"""
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    lang = config.lang.rscript
    envs = {
        "infmt": "matrix",  # or rds
        "inunit": None,
        "outunit": None,
        "refexon": config.ref.refexon,
        "meanfl": None,
        "inlog2p": False,
        "outlog2p": False,
    }
    script = "file://../scripts/rnaseq/UnitConversion.R"


class Simulation(Proc):
    """Simulate RNA-seq data using ESCO package

    Input:
        ngenes: Number of genes to simulate
        nsamples: Number of samples to simulate
            If you want to force the process to re-simulate for the same
            `ngenes` and `nsamples`, you can set a different value for `envs.seed`.
            Note that the samples will be shown as cells in the output (since
            the simulation is designed for single-cell RNA-seq data).

    Output:
        outdir: Output directory containing the simulated data

    Envs:
        ncores (type=int): Number of cores to use.
        type (choice): Which type of heterogenounity to use.
            - single: produces a single population.
            - group: produces distinct groups.
            - tree: produces distinct groups but admits a tree structure.
            - traj: produces distinct groups but admits a smooth trajectory structure.
        seed (type=int): Random seed.
            If not set, seed will not be set.
        args (type=json): Additional arguments to pass to the simulation function.
            See <https://rdrr.io/github/JINJINT/ESCO/man/escoParams.html>.
    """
    input = "ngenes:var, nsamples:var"
    output = "outdir:dir:{{in.ngenes}}x{{in.nsamples}}.sim"
    lang = config.lang.rscript
    envs = {
        "ncores": config.misc.ncores,
        "type": "single",
        "args": {"dropout-type": "none"},
        "seed": None,
    }
    script = "file://../scripts/rnaseq/Simulation.R"
