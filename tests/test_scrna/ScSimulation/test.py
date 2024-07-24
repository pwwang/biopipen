from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import ScSimulation as ScSimulation_
from biopipen.core.testing import get_pipeline
from datar.tibble import tibble


class ScSimulation1(ScSimulation_):
    input_data = [1234]


class ScSimulation2(ScSimulation_):
    input_data = ["Sample1"]
    envs = {"params": {"group-prob": [0.3, 0.7]}, "method": "groups"}


def pipeline():
    return get_pipeline(__file__).set_starts(ScSimulation1, ScSimulation2)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
