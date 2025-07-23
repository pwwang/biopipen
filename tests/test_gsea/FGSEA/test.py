from biopipen.core.config import config
from biopipen.core.proc import Proc
from biopipen.ns.web import Download
from biopipen.ns.misc import Str2File
from biopipen.ns.gsea import FGSEA as FGSEA_
from biopipen.core.testing import get_pipeline


class Preprocessing(Proc):
    requires = Download
    input = "infile:file"
    input_data = lambda ch: [ch.outfile[0]]
    output = "outfile:file:{{in.infile | stem}}.txt"
    lang = config.lang.rscript
    script = """
        infile <- {{in.infile | quote}}
        outfile <- {{out.outfile | quote}}
        data <- read.table(infile, header = TRUE, row.names = NULL, sep = "\\t", check.names = FALSE)
        colnames(data)[1] <- "Gene"
        data <- data[!duplicated(data$Gene), ]
        rownames(data) <- data$Gene
        data$Gene <- NULL
        write.table(data, outfile, sep = "\\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    """


FGSEA = Proc.from_proc(
    FGSEA_,
    requires=[Preprocessing, Str2File],
    input_data=lambda ch1, ch2: [(ch1.outfile[0], ch2.outfile[0])],
    envs={
        "clscol": "Group",
        "case": "MMP9",
        "control": "CTRL",
        "gmtfile": "KEGG_2021_Human",
    },
)


def pipeline():
    return (
        get_pipeline(__file__)
        .set_starts(
            Download,
            Str2File,
        )
        .set_data(
            [
                "https://www.ncbi.nlm.nih.gov/geo/download/"
                "?acc=GSE179367"
                "&format=file"
                "&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz",
            ],
            [
                (
                    (
                        "Sample\tGroup\n"
                        "shMMP9_1\tMMP9\n"
                        "shMMP9_2\tMMP9\n"
                        "shMMP9_3\tMMP9\n"
                        "shControl_1\tCTRL\n"
                        "shControl_2\tCTRL\n"
                        "shControl_3\tCTRL\n"
                    ),
                    "samples.txt",
                )
            ],
        )
    )


def testing(pipen):
    # assert pipen._succeeded
    outfile = pipen.procs[-1].workdir.joinpath(
        "0",
        "output",
        "acc.GSE179367.format.file.file.GSE179367_gene_count.real.txt.fgsea",
    )
    assert outfile.is_dir()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
