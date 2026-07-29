from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import (
    ScRepLoading as ScRepLoading_,
    ScRepCombiningExpression as ScRepCombiningExpression_,
    ClonalStats as ClonalStats_,
)
from biopipen.core.testing import get_pipeline


class ContigPreparation(Proc):
    """Prepare the contig files for ScRepLoading"""

    input = "seed:var"
    input_data = [8525]
    output = "metafile:file:metafile.tsv, contigdir:dir:contigs"
    lang = config.lang.rscript
    script = """
        library(scRepertoire)
        metafile <- {{out.metafile | quote}}
        contigdir <- {{out.contigdir | quote}}

        data(contig_list)
        names(contig_list) <- c("P17B", "P17L", "P18B", "P18L",
                                "P19B", "P19L", "P20B", "P20L")
        lapply(names(contig_list), function(name) {
            contigs <- contig_list[[name]]
            ctgdir <- file.path(contigdir, name)
            dir.create(ctgdir, showWarnings = FALSE)
            write.csv(contigs, row.names = FALSE,
                file.path(ctgdir, "filtered_contig_annotations.csv"))
        })

        write.table(
            data.frame(Sample = names(contig_list),
                Type = rep(c("B", "L"), 4),
                TCRData = file.path(
                    contigdir, names(contig_list), "filtered_contig_annotations.csv")),
            metafile, sep = "\\t", row.names = FALSE, quote = FALSE
        )
    """


class ScRepLoading(ScRepLoading_):
    requires = ContigPreparation


class ScExpression(Proc):
    """Load the expression data for ScRepCombiningExpression"""

    input = "seed:var"
    input_data = [8525]
    output = "exprfile:file:exprfile.rds"
    lang = config.lang.rscript
    script = """
        library(Seurat)
        library(scRepertoire)
        exprfile <- {{out.exprfile | quote}}

        data(scRep_example)
        scRep_example$Sample <- scRep_example$orig.ident
        scRep_example <- UpdateSeuratObject(scRep_example)
        saveRDS(scRep_example, exprfile)
    """


class ScRepCombiningExpression(ScRepCombiningExpression_):
    requires = ScRepLoading, ScExpression


class ClonalStats(ClonalStats_):
    requires = ScRepCombiningExpression
    envs_depth = 2
    envs = {
        "cases": {
            "CDR3 Length by cloneSize": {
                "viz_type": "length",
                "group_by": "cloneSize",
                "plot_type": "box",
            },
            "Clonal Volume": {
                "viz_type": "volume",
                "save_data": True,
            },
            "Clonal Abundance": {"viz_type": "abundance"},
        }
    }


def pipeline():
    return get_pipeline(__file__, enable_report=False).set_starts(
        ContigPreparation, ScExpression
    )


def testing(pipen): ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
