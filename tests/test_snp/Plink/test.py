from pathlib import Path

from biopipen.ns.misc import File2Proc
from biopipen.ns.snp import (
    PlinkFromVcf as PlinkFromVcf_,
    Plink2GTMat as Plink2GTMat_,
    PlinkIBD as PlinkIBD_,
    PlinkHWE as PlinkHWE_,
    PlinkHet as PlinkHet_,
    PlinkCallRate as PlinkCallRate_,
    PlinkFilter as PlinkFilter_,
    PlinkFreq as PlinkFreq_,
    PlinkUpdateName as PlinkUpdateName_,
)
from biopipen.core.testing import get_pipeline


class VCFFile(File2Proc):
    ...


class NameFile1(File2Proc):
    ...


class NameFile2(File2Proc):
    ...


class PlinkFromVcf1(PlinkFromVcf_):
    requires = VCFFile


class PlinkFromVcf2(PlinkFromVcf_):
    requires = VCFFile
    envs = {"set_missing_var_ids": None}


class PlinkUpdateName1(PlinkUpdateName_):
    requires = PlinkFromVcf1, NameFile1


class PlinkUpdateName2(PlinkUpdateName_):
    requires = PlinkFromVcf1, NameFile2


class PlinkFilter(PlinkFilter_):
    requires = PlinkFromVcf1
    envs = {"samples": ["A", "B"], "keep": True}


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


class PlinkFreqNoFilter(PlinkFreq_):
    requires = PlinkFromVcf1
    envs = {"cutoff": 0.2}


class PlinkFreqGtFilter(PlinkFreq_):
    requires = PlinkFromVcf1
    envs = {"cutoff": 0.2, "filter": "gt"}


class PlinkFreqCountsNoFilter(PlinkFreq_):
    requires = PlinkFromVcf1
    envs = {"cutoff": 3, "modifier": "counts"}


class PlinkFreqCountsLtFilter(PlinkFreq_):
    requires = PlinkFromVcf1
    envs = {"cutoff": 3, "filter": "lt", "modifier": "counts"}


# class PlinkFreqCCNoFilter(PlinkFreq_):
#     requires = PlinkFromVcf1
#     envs = {"cutoff": 0.2, "modifier": "cc"}


# class PlinkFreqCCLeFilter(PlinkFreq_):
#     requires = PlinkFromVcf1
#     envs = {
#         "cutoff": {"NCHROBS_A": 0, "NCHROBS_U": 0},
#         "filter": "lt",
#         "modifier": "cc",
#     }


class PlinkFreqXLeFilter(PlinkFreq_):
    requires = PlinkFromVcf1
    envs = {
        "cutoff": {"HET_REF_ALT1_CT": 2, "HOM_ALT1_CT": 2},
        "filter": "lt",
        "modifier": "x",
    }


class PlinkCallRate1(PlinkCallRate_):
    requires = PlinkFromVcf1


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["+report"])
        get_pipeline(__file__)
        .set_starts(VCFFile, NameFile1, NameFile2)
        .set_data(
            [Path(__file__).parent / "data" / "sample.vcf.gz"],
            [Path(__file__).parent / "data" / "name_to_rs.txt"],
            [Path(__file__).parent / "data" / "name_to_rs.vcf.gz"],
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
