from datar.tibble import tibble
from pipen import Pipen, Proc
from biopipen.core.config import config
from biopipen.namespaces.scrna import SeuratSubset as SeuratSubset_


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


class SeuratSubset(SeuratSubset_):
    requires = PrepareSeurat
    input_data = lambda ch: tibble(
        srtobj=ch,
        filters=[
            {"Group": {"mutaters": {"Group": "groups"}, "groupby": "Group"}}
        ],
    )


pipen = (
    Pipen("TestScrnaSeuratSubset")
    .set_start(PrepareSeurat)
)

if __name__ == "__main__":
    pipen.run()
