from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import (
    ScRepLoading as ScRepLoading_,
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


class ClonalStats(ClonalStats_):
    requires = ScRepLoading
    envs_depth = 2
    envs = {
        "cases": {
            "CDR3 Length by Type": {
                "viz_type": "length",
                "group_by": "Type",
                "plot_type": "box",
            }
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(ContigPreparation)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
