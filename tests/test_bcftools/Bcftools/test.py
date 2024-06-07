from pathlib import Path

from biopipen.ns.bcftools import (
    BcftoolsAnnotate as BcftoolsAnnotate_,
    BcftoolsFilter as BcftoolsFilter_,
    BcftoolsView as BcftoolsView_,
    BcftoolsSort as BcftoolsSort_,
)
from biopipen.core.testing import get_pipeline

vcffile = Path(__file__).parent.parent.parent.joinpath(
    "test_snp", "Plink", "data", "sample.vcf.gz"
)
annofile = Path(__file__).parent.parent.joinpath("Bcftools", "data", "anno.vcf")
regfile = Path(__file__).parent.parent.joinpath("Bcftools", "data", "regions.bed")


class BcftoolsAnnotateRemove(BcftoolsAnnotate_):
    input_data = [vcffile]
    envs = {"remove": ["INFO/VCB", "FORMAT/DV"]}


class BcftoolsSort(BcftoolsSort_):
    input_data = [annofile]
    envs = {"gz": False, "index": False}


class BcftoolsAnnotateAdd(BcftoolsAnnotate_):
    requires = [BcftoolsAnnotateRemove, BcftoolsSort]
    envs = {"columns": ["INFO/ADP", "INFO/AUGT"]}


class BcftoolsFilter(BcftoolsFilter_):
    input_data = [vcffile]
    envs = {
        "includes": "MIN(FMT/DP)>10",
        "excludes": "INFO/DP>500",
    }


class BcftoolsView(BcftoolsView_):
    input_data = [vcffile]
    envs = {
        "regions_file": regfile,
    }


def pipeline():
    return (
        # get_pipeline(__file__, plugins=["+report"])
        get_pipeline(__file__)
        .set_starts(
            BcftoolsAnnotateRemove,
            BcftoolsFilter,
            BcftoolsView,
            BcftoolsSort,
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
