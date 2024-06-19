from pathlib import Path

from biopipen.ns.stats import (
    MetaPvalue as MetaPvalue_,
    MetaPvalue1 as MetaPvalue1_,
)
from biopipen.core.testing import get_pipeline

infile1 = Path(__file__).parent / "data" / "pvalues1.txt"
infile2 = Path(__file__).parent / "data" / "pvalues2.txt"
infile3 = Path(__file__).parent / "data" / "pvalues3.txt"
infile = Path(__file__).parent / "data" / "pvalues-single.txt"


class MetaPvalue(MetaPvalue_):
    envs = {
        "id_cols": "ID1",
        "pval_cols": "pval",
        "padj": "BH",
    }


class MetaPvalue_1(MetaPvalue_):
    envs = {
        "id_cols": "ID",
        "id_exprs": ["ID1", "ID3"],
        "pval_cols": ["pval", "pval3"],
        "padj": "none",
    }


class MetaPvalue1(MetaPvalue1_):
    envs = {
        "id_cols": "ID1,ID2",
        "pval_col": "pval",
        "padj": "fdr",
    }


def pipeline():
    return (
        get_pipeline(__file__)
        .set_starts(MetaPvalue, MetaPvalue_1, MetaPvalue1)
        .set_data(
            [[infile1, infile2]],
            [[infile1, infile3]],
            [infile],
        )
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
