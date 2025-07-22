from biopipen.ns.rnaseq import Simulation as Simulation_
from biopipen.core.testing import get_pipeline


class SimulationRUVcorr(Simulation_):
    input_data = [(100, 10)]


class SimulationESCO(Simulation_):
    input_data = [(100, 10)]
    envs = {"tool": "ESCO"}


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__)
        .set_starts(SimulationRUVcorr, SimulationESCO)
    )


def testing(pipen):
    # assert pipen._succeeded
    pass


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
