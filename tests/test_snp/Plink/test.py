from pathlib import Path

from biopipen.ns.misc import File2Proc
from biopipen.ns.snp import (
    PlinkFromVcf as PlinkFromVcf_,
    Plink2GTMat as Plink2GTMat_,
    PlinkIBD as PlinkIBD_,
    PlinkHWE as PlinkHWE_,
    PlinkHet as PlinkHet_,
    PlinkCallRate as PlinkCallRate_,
)
from biopipen.core.testing import get_pipeline


class PlinkFromVcf1(PlinkFromVcf_):
    requires = File2Proc


class PlinkFromVcf2(PlinkFromVcf_):
    requires = File2Proc
    envs = {"set_missing_var_ids": None}


class Plink2GTMat1(Plink2GTMat_):
    requires = PlinkFromVcf1


class Plink2GTMat2(Plink2GTMat_):
    requires = PlinkFromVcf2
    envs = {"transpose": True}


class PlinkIBD1(PlinkIBD_):
    requires = PlinkFromVcf1


class PlinkHWE1(PlinkHWE_):
    requires = PlinkFromVcf1
    envs = {"filter": True, "cutoff": 0.6}


class PlinkHet1(PlinkHet_):
    requires = PlinkFromVcf1
    envs = {"filter": True, "cutoff": 1}


class PlinkCallRate1(PlinkCallRate_):
    requires = PlinkFromVcf1


def pipeline():
    return (
        get_pipeline(__file__, plugins=["+report"])
        .set_starts(File2Proc)
        .set_data([Path(__file__).parent / "data" / "sample.vcf.gz"])
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
