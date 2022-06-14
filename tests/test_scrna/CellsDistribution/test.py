from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import CellsDistribution
from biopipen.core.testing import get_pipeline


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
        set.seed(8525)
        pbmc_small@meta.data$response = sample(
            1:3,
            nrow(pbmc_small@meta.data),
            replace=TRUE
        )
        pbmc_small@meta.data$Clone = paste0(
            "C",
            sample(1:10, nrow(pbmc_small@meta.data), replace=TRUE)
        )
        pbmc_small$Sample = pbmc_small$letter.idents
        pbmc_small$seurat_clusters = pbmc_small$groups
        saveRDS(pbmc_small, {{out.outfile | quote}})
    """


class CellsDistribution(CellsDistribution):
    requires = PrepareSeurat
    envs = {
        "cases": {
            "cases": {
                "Case1": {
                    "group": {
                        "by": "Responder",
                        "mutaters": {
                            "Responder": """
                                case_when(
                                    response == 2 ~ "Responder",
                                    response == 3 ~ "NonResponder",
                                    TRUE ~ "CONTROL"
                                )
                            """
                        },
                        "order": ["CONTROL", "Responder", "NonResponder"],
                    },
                    "cells": {
                        "by": "Clone",
                        "n": 5,
                        "orderby": "desc(.CloneSize)"
                    },
                }
            }
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareSeurat)


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath(
            "0",
            "output",
            "pbmc_small.cells_distribution",
        )
    )
    assert outfile.is_dir()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
