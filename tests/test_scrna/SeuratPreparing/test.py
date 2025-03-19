from pathlib import Path

from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import (
    SeuratTo10X as SeuratTo10X_,
    SeuratPreparing as SeuratPreparing_,
)
from biopipen.core.testing import get_pipeline, _find_testing_index, TESTING_DIR


class PrepareSeurat(Proc):
    """Prepare the data"""
    input = "name"
    input_data = ["pbmc3k"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        options(timeout=600)
        library(SeuratData)
        InstallData("pbmc3k")
        pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
        pbmc3k$Sample <- paste0("S", sample(1:2, nrow(pbmc3k), replace = TRUE))
        saveRDS(pbmc3k, {{out.outfile | quote}})
    """


class SeuratTo10X(SeuratTo10X_):
    requires = PrepareSeurat
    envs = {"split_by": "Sample"}


def gen_input(ch):
    metafile = Path(TESTING_DIR % {"index": _find_testing_index(False)}).joinpath(
        "SeuratPreparing-test.txt"
    )
    if not metafile.is_file():
        with metafile.open("w") as f:
            f.write("Sample\tRNAData\n")
            for s in Path(ch.iloc[0, 0]).glob("*"):
                f.write(f"{s.name}\t{s}\n")

    return [metafile]


class SeuratPreparing(SeuratPreparing_):
    requires = SeuratTo10X
    input_data = gen_input
    envs = {"cell_qc": "runif(n()) < 0.5", "IntegrateLayers": {"method": "rpca"}}


class SeuratPreparing2(SeuratPreparing_):
    requires = SeuratTo10X
    input_data = gen_input
    envs = {
        "cell_qc": "runif(n()) < 0.5",
        "cell_qc_per_sample": True,
        "doublet_detector": "DoubletFinder",
        "DoubletFinder": {"PCs": 3},
    }


class SeuratPreparing3(SeuratPreparing_):
    requires = SeuratTo10X
    input_data = gen_input
    envs = {
        "cell_qc": "runif(n()) < 0.5",
        "cell_qc_per_sample": True,
        "doublet_detector": "scDblFinder",
    }


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_starts(PrepareSeurat)
    )


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
