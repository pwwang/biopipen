from pathlib import Path

from biopipen.ns.bed import Bed2Vcf
from biopipen.core.testing import get_pipeline

class Bed2Vcf(Bed2Vcf):

    envs = {
        "infos": [
            {
                "ID": "SVType",
                "Number": "1",
                "Type": "String",
                "Description": "Type of structural variant"
            },
        ],
        "formats": [
            {
                "ID": "Depth",
                "Number": "1",
                "Type": "Interger",
                "Description": "Depth of coverage"
            },
        ],
        "converters": {
            "ID": "lambda items: items[3]",
            "ALT": "lambda items: f'<{items[4]}>'",
            "SVType": "lambda items: items[4]",
            "Depth": "lambda items: int(items[5])",
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(Bed2Vcf).set_data(
        [Path(__file__).parent / "data" / "in.bed"]
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "in.vcf.gz")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
