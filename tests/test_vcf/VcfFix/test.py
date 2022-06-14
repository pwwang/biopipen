from pathlib import Path

from biopipen.ns.vcf import VcfFix
from biopipen.core.testing import get_pipeline


class VcfFix(VcfFix):

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


INFILE = Path(__file__).parent / "data" / "tofix.vcf"


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(VcfFix)
        .set_data([str(INFILE)])
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", INFILE.name)
    )
    assert outfile.is_file()

if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
