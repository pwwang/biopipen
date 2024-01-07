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
        pbmc_small = FindClusters(pbmc_small)
        saveRDS(pbmc_small, {{out.outfile | quote}})
    """


class ScFGSEA(ScFGSEA):
    requires = PrepareSeurat
    envs = {
        # "each": "seurat_clusters",
        "ident-1": "g1",
        "ident-2": "g2",
        "group-by": "groups",
        # "gmtfile": Path(__file__).parent.parent.parent.joinpath(
        #     "data/reference/KEGG_metabolism.gmt"
        # ),
        "subset": "seurat_clusters == '0'",
        "gmtfile": (
            "https://download.baderlab.org/EM_Genesets/January_01_2024/"
            "Human/symbol/Pathways/Human_Hallmark_MSigdb_January_01_2024_symbol.gmt"
        ),
    }


def pipeline():
    return (
        get_pipeline(__file__, plugins=["no:report"])
        .set_start(PrepareSeurat)
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc_small.fgsea"
        )
    )
    assert outfile.is_dir()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
