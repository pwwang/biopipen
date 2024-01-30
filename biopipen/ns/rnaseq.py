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
    """Simulate RNA-seq data using ESCO/RUVcorr package

    Input:
        ngenes: Number of genes to simulate
        nsamples: Number of samples to simulate
            If you want to force the process to re-simulate for the same
            `ngenes` and `nsamples`, you can set a different value for `envs.seed`.
            Note that the samples will be shown as cells in the output (since
            the simulation is designed for single-cell RNA-seq data).

    Output:
        outdir: Output directory containing the simulated data
            `sim.rds` and `True.rds` will be generated.
            For `ESCO`, `sim.rds` contains the simulated data in a
            `SingleCellExperiment` object, and `True.rds` contains the matrix of true
            counts with columns
            representing samples and rows representing genes.
            For `RUVcorr`, `sim.rds` contains the simulated data in list with
            `Truth`, A matrix containing the values of Xβ; `Y` A matrix containing the
            values in `Y`; `Noise` A matrix containing the values in `Wα`; `Sigma`
            A matrix containing the true gene-gene correlations, as defined by Xβ; and
            `Info` A matrix containing some of the general information about the
            simulation. For the matrices, rows represent genes and columns represent
            samples.

    Envs:
        tool (choice): Which tool to use for simulation.
            - ESCO: uses the [ESCO](https://github.com/JINJINT/ESCO) package.
            - RUVcorr: uses the [RUVcorr](https://rdrr.io/bioc/RUVcorr/) package.
        ncores (type=int): Number of cores to use.
        seed (type=int): Random seed.
            If not set, seed will not be set.
        esco_args (ns): Additional arguments to pass to the simulation function.
            - type (choice): Which type of heterogenounity to use.
                - single: produces a single population.
                - group: produces distinct groups.
                - tree: produces distinct groups but admits a tree structure.
                - traj: produces distinct groups but admits a smooth trajectory
                    structure.
            - <more>: See <https://rdrr.io/github/JINJINT/ESCO/man/escoParams.html>.
        ruvcorr_args (ns): Additional arguments to pass to the simulation
            function.
            - <more>: See <https://rdrr.io/bioc/RUVcorr/man/simulateGEdata.html>.
    """
    input = "ngenes:var, nsamples:var"
    output = "outdir:dir:{{in.ngenes}}x{{in.nsamples}}.sim"
    lang = config.lang.rscript
    envs = {
        "tool": "RUVcorr",
        "ncores": config.misc.ncores,
        "type": "single",
        "esco_args": {"dropout-type": "none"},
        "ruvcorr_args": {},
        "seed": None,
    }
    script = "file://../scripts/rnaseq/Simulation.R"
