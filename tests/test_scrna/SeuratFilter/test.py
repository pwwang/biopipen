from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import SeuratFilter
from biopipen.core.testing import get_pipeline
from datar.tibble import tibble


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


class SeuratFilter(SeuratFilter):
    requires = PrepareSeurat
    input_data = lambda ch: tibble(
        srtobj=ch,
        filters=[{"mutaters": {"GROUP": "groups"}, "filter": "GROUP == 'g1'"}],
    )


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareSeurat)


def testing(pipen):
    assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc_small.filtered.RDS",
        )
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
