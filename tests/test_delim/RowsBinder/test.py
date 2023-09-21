from pathlib import Path

from biopipen.ns.delim import RowsBinder as RowsBinder_
from biopipen.core.testing import get_pipeline


class RowsBinder(RowsBinder_):
    ...


class RowsBinderWithFilenames(RowsBinder_):
    envs = {"filenames": ["A", "B"]}


class RowsBinderWithFilenameStrings(RowsBinder_):
    envs = {"filenames": "A, B"}


class RowsBinderWithFilenameFuns(RowsBinder_):
    envs = {"filenames": "function(x) tools::file_path_sans_ext(basename(x))"}


def pipeline():
    return (
        get_pipeline(__file__, plugins=["no:report"])
        .set_starts(
            RowsBinder,
            RowsBinderWithFilenames,
            RowsBinderWithFilenameStrings,
            RowsBinderWithFilenameFuns,
        )
        .set_data(
            *[
                [
                    [
                        Path(__file__).parent / "data" / "A.txt",
                        Path(__file__).parent / "data" / "B.txt",
                    ]
                ]
            ] * 4
        )
    )


def testing(pipen):
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "A_rbound.txt")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
