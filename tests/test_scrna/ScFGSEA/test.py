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


class ScFGSEASingle(ScFGSEA):
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
            "https://raw.githubusercontent.com/pwwang/immunopipe-example/"
            "master/data/KEGG_metabolism.short.gmt"
        ),
    }
    order = 9


class ScFGSEAEach(ScFGSEA):
    requires = PrepareSeurat
    envs = {
        # Test if seurat_clusters == 1 gets ignored
        "mutaters": {"groups2": "ifelse(seurat_clusters == '1', NA, groups)"},
        "ident-1": "g1",
        "ident-2": "g2",
        "group-by": "groups2",
        # "gmtfile": Path(__file__).parent.parent.parent.joinpath(
        #     "data/reference/KEGG_metabolism.gmt"
        # ),
        "each": "seurat_clusters",
        "gmtfile": (
            "https://raw.githubusercontent.com/pwwang/immunopipe-example/"
            "master/data/KEGG_metabolism.short.gmt"
        ),
    }


def pipeline():
    return get_pipeline(__file__).set_start(PrepareSeurat)


def testing(pipen):
    # assert pipen._succeeded
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
    assert pipen.run()
    testing(pipen)
