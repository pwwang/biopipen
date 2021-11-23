from pipen import Pipen, Proc
from biopipen.namespaces.misc import Str2File
from biopipen.namespaces.gene import GeneNameConversion

GeneNameConversion1 = Proc.from_proc(
    GeneNameConversion,
    requires=Str2File,
    envs={"infmt": "reporter", "genecol": 1, "notfound": "use-query"},
)

GeneNameConversion2 = Proc.from_proc(
    GeneNameConversion,
    requires=Str2File,
    envs={
        "infmt": "reporter",
        "outfmt": "ensembl.gene,symbol",
        "notfound": "skip",
        "genecol": 1,
    },
)

GeneNameConversion3 = Proc.from_proc(
    GeneNameConversion,
    requires=Str2File,
    envs={
        "infmt": "reporter",
        "outfmt": "ensembl.gene,symbol",
        "output": "replace",
        "notfound": "skip",
        "genecol": 1,
    },
)

# SHould raise error
# GeneNameConversion4 = Proc.from_proc(
#     GeneNameConversion,
#     requires=Str2File,
#     envs={"infmt": "reporter", "genecol": 1, "notfound": "error"},
# )

pipen = (
    Pipen("TestGeneGeneNameConversion")
    .set_start(Str2File)
    .set_data(
        [
            (
                (
                    "Meta\tId\tsymbol\n"
                    "a\t1255_g_at\t1\n"
                    "b\t1294_at\t2\n"
                    "c\t1316_at\t3\n"
                    "d\tnonexist\t4\n"
                    "e\t1320_at\t5\n"
                ),
                "data.txt",
            )
        ]
    )
)

if __name__ == "__main__":
    pipen.run()
