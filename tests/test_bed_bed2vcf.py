from pathlib import Path

from pipen import Pipen
from biopipen.namespaces.bed import Bed2Vcf as Bed2Vcf_

class Bed2Vcf(Bed2Vcf_):

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


pipen = Pipen().set_starts(Bed2Vcf).set_data(
    [Path(__file__).parent / "data" / "bed" / "Bed2Vcf" / "in.bed"]
)

if __name__ == "__main__":
    pipen.run()
