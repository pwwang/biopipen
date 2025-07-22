import os
from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import CDR3AAPhyschem as CDR3AAPhyschem_, TESSA as TESSA_
from biopipen.core.testing import get_pipeline

conda_exe = os.environ.get("CONDA_EXE", "")
if conda_exe:
    # Make sure we have scvelo 0.3.3 and numpy<2
    python = os.path.join(
        os.path.dirname(os.path.dirname(conda_exe)), "bin", "python"
    )
else:
    python = "python"


class DataPreparation(Proc):
    """Prepare the data files for CDR3AAPhyschem"""
    input = "seed:var"
    input_data = [8525]
    output = "datafile:file:data.qs"
    lang = config.lang.rscript
    script = """
        library(scRepertoire)
        library(biopipen.utils)
        set.seed({{in.seed | int}})

        # Getting the combined contigs
        combined <- combineTCR(contig_list,
                                samples = c("P17B", "P17L", "P18B", "P18L",
                                            "P19B","P19L", "P20B", "P20L"))

        # Getting a sample of a Seurat object
        scRep_example <- get(data("scRep_example"))

        # Using combineExpresion()
        scRep_example <- combineExpression(combined, scRep_example)

        # Save the data
        biopipen.utils::save_obj(scRep_example, {{out.datafile | quote}})
    """


class CDR3AAPhyschem(CDR3AAPhyschem_):
    requires = DataPreparation
    envs = {
        "group": "seurat_clusters",
        "comparison": {
            "Treg": ["3", "4", "5"],
            "Tconv": ["6", "7", "8"],
        }
    }


class TESSA(TESSA_):
    requires = DataPreparation
    envs = {"max_iter": 100, "python": python}


def pipeline():
    return get_pipeline(__file__).set_starts(DataPreparation)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
