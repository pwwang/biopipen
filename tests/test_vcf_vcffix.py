from pathlib import Path

from pipen import Pipen
from biopipen.ns.vcf import VcfFix as VcfFix_


class VcfFix(VcfFix_):

    envs = {
        "fixes": [
            {
                "kind": "contig",
                "append": True,
                "fix": "lambda obj: HeaderContig(ID='chr1', length=249250621)",
            },
            {
                "kind": "format",
                "id": "Depth",
                "fix": "lambda obj: setattr(obj, 'Type', 'String')",
            },
            {
                "kind": "fields",
                "fix": "lambda items: items.append(instem)",
            },
            {
                "kind": "variant",
                "append": True,
                "fix": """lambda items: Variant.from_strs(
                    'chr1',
                    50,
                    'S3',
                    'A',
                    '<DEL>',
                    '.',
                    'PASS',
                    'SVType=DEL;END=60',
                    'GT:Depth',
                    '0/1:3',
                )""",
            },
        ],
    }


pipen = Pipen().set_starts(VcfFix).set_data(
    [Path(__file__).parent / "data" / "vcf" / "VcfFix" / "tofix.vcf"]
)

if __name__ == "__main__":
    pipen.run()
