import os
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import ScVelo as ScVelo_
from biopipen.core.testing import get_pipeline

conda_prefix = os.environ.get("CONDA_PREFIX", "")
if conda_prefix:
    # Make sure we have scvelo 0.3.3 and numpy<2
    python = os.path.join(conda_prefix, "envs", "bio39", "bin", "python")
else:
    python = "python"


class PrepareAnnData(Proc):
    lang = python
    input = "var"
    input_data = ["anndata"]
    output = "outfile:file:{{in.var}}.h5ad"
    script = """
        import scanpy as sc
        import scvelo as scv

        adata = scv.datasets.pancreas()
        adata.write_h5ad({{out.outfile | quote}})
    """


class PrepareSeurat(Proc):
    requires = PrepareAnnData
    lang = config.lang.rscript
    input = "annfile:file"
    output = "outfile:file:{{in.annfile | stem}}.qs"
    script = """
        biopipen.utils::ConvertAnnDataToSeurat({{in.annfile | r}}, {{out.outfile | r}})
    """


class ScVeloAnnData(ScVelo_):
    requires = PrepareAnnData
    lang = python
    envs = {"group_by": "clusters", "ncores": 16}


class ScVeloSeurat(ScVelo_):
    requires = PrepareSeurat
    lang = python
    envs = {"group_by": "clusters", "ncores": 16}


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareAnnData)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
