from pipen import Pipen, Proc
from biopipen.namespaces.web import Download
from biopipen.namespaces.misc import Str2File
from biopipen.namespaces.gsea import FGSEA
from datar.all import flatten, tibble, select

FGSEA = Proc.from_proc(
    FGSEA,
    requires=[Download, Str2File],
    input_data=lambda ch1, ch2: tibble(*flatten(ch1), ch2) >> select(1, 3, 2),
    envs={"clscol": "Group", "classes": ["MMP9", "CTRL"]}
)

pipen = (
    Pipen("TestGseaFGSEA")
    .set_start(Download, Str2File)
    .set_data(
        [
            "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE179367&format=file&file=GSE179367%5Fgene%5Fcount%2Ereal%2Etxt%2Egz",
            "https://www.genepattern.org/tutorial/linkedFiles/export_gnf.GENE_SYMBOL.gmt",
        ],
        [(
            (
                "Sample\tGroup\n"
                "shMMP9_1\tMMP9\n"
                "shMMP9_2\tMMP9\n"
                "shMMP9_3\tMMP9\n"
                "shControl_1\tCTRL\n"
                "shControl_2\tCTRL\n"
                "shControl_3\tCTRL\n"
            ),
            "samples.txt"
        )]
    )
)

if __name__ == "__main__":
    pipen.run()
