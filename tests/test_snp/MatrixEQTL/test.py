from pathlib import Path
from biopipen.ns.snp import MatrixEQTL as MatrixEQTL_
from biopipen.core.testing import get_pipeline

datadir = Path(__file__).parent / "data"


class MatrixEQTL(MatrixEQTL_):
    input_data = [
        (datadir / "SNP.txt", datadir / "GE.txt", datadir / "Covariates.txt")
    ]


class MatrixEQTLNoCov(MatrixEQTL_):
    input_data = [
        (datadir / "SNP.txt", datadir / "GE.txt", None)
    ]


class MatrixEQTLCis(MatrixEQTL_):
    input_data = [
        (datadir / "SNP.txt", datadir / "GE.txt", datadir / "Covariates.txt")
    ]
    envs = {"snppos": datadir / "snpsloc.txt", "genepos": datadir / "geneloc.txt"}


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__)
        .set_starts(MatrixEQTL, MatrixEQTLNoCov, MatrixEQTLCis)
    )


def testing(pipen):
    # assert pipen._succeeded
    pass


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
