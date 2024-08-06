from biopipen.ns.misc import Str2File
from biopipen.ns.gene import GeneNameConversion
from biopipen.core.testing import get_pipeline


class InputFile(Str2File):
    input_data = [
        (
            "Meta\tEnsg\n"
            # "a\tENSG00000230373\n"
            "b\tENSG00000236269\n"
            "c\tENSG00000227232.5\n"
            "c1\tENSG00000227232.5\n"
            "d\tnonexist\n"
        )
    ]


class GeneNameConversion1(GeneNameConversion):
    requires = InputFile
    envs = {
        "infmt": "ensg",
        "outfmt": "symbol",
        "notfound": "use-query",
        "output": "append",
        "genecol": 2,
    }


class GeneNameConversion2(GeneNameConversion):
    requires = InputFile
    envs = {
        "infmt": "ensg",
        "outfmt": "symbol",
        "notfound": "na",
        "output": "replace",
        "genecol": 2,
    }


class GeneNameConversion3(GeneNameConversion):
    requires = InputFile
    envs = {
        "infmt": "ensg",
        "outfmt": "symbol",
        "notfound": "na",
        "output": "with-query",
        "dup": "combine",
        "genecol": 2,
    }


def pipeline():
    return get_pipeline(__file__).set_starts(InputFile)


def _assert_output(seq1: str, seq2: str, proc: str):
    seq1 = seq1.splitlines()
    seq2 = seq2.splitlines()
    assert len(seq1) == len(seq2), f"{proc} output length mismatch"
    for i, (s1, s2) in enumerate(zip(seq1, seq2)):
        assert s1 == s2, f"{proc} output line {i + 1} mismatch: {s1!r} != {s2!r}"


def _testing_proc(proc: str, exp: str):
    out = pipeline().outdir.joinpath(proc, "unnamed.txt")
    assert out.is_file(), out
    out = out.read_text()
    _assert_output(out, exp, proc)


def testing(pipen):
    exp1 = (
        "Meta\tEnsg\tsymbol\n"
        # "a\tENSG00000230373\tGOLGA6L17P\n"
        "b\tENSG00000236269\tENSG00000236269\n"
        "c\tENSG00000227232.5\tWASH7P\n"
        "c1\tENSG00000227232.5\tWASH7P\n"
        "d\tnonexist\tnonexist\n"
    )
    _testing_proc("GeneNameConversion1", exp1)

    exp2 = (
        "Meta\tEnsg\n"
        # "a\tGOLGA6L17P\n"
        "b\tNA\n"
        "c\tWASH7P\n"
        "c1\tWASH7P\n"
        "d\tNA\n"
    )
    _testing_proc("GeneNameConversion2", exp2)

    exp3 = (
        "query\tsymbol\n"
        # "ENSG00000230373\tGOLGA6L17P;GOLGA6L3P\n"
        "ENSG00000236269\tNA\n"
        "ENSG00000227232.5\tWASH7P\n"
        "ENSG00000227232.5\tWASH7P\n"
        "nonexist\tNA\n"
    )
    _testing_proc("GeneNameConversion3", exp3)


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
