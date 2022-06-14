from pathlib import Path

from pipen import Proc
from biopipen.core.config import config
from biopipen.core.testing import get_pipeline
from biopipen.ns.scrna import ScFGSEA


class PrepareSeurat(Proc):
    """Prepare the data

    Requires:
        - name: Seurat
          check: |
            {{proc.lang}} <(echo "library(Seurat)")
    """

    input = "name"
    input_data = ["pbmc_small"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        library(Seurat)
        saveRDS(pbmc_small, {{out.outfile | quote}})
    """


class ScFGSEA(ScFGSEA):
    requires = PrepareSeurat
    envs = {
        "cases": {
            "name": "GSEA analysis",
            "cases": {
                "Male_vs_Female": {
                    "ident.1": "g1",
                    "ident.2": "g2",
                    "group.by": "Group",
                    "mutaters": {"Group": "groups"},
                }
            }
        },
        "gmtfile": Path(__file__).parent.parent.parent.joinpath(
            "data/reference/KEGG_metabolism.gmt"
        ),
    }


def pipeline():
    return (
        get_pipeline(__file__)
        .set_start(PrepareSeurat)
    )


def testing(pipen):
    outfile = (
        pipen.outdir.joinpath(
            "REPORTS",
            "index.html",
        )
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
