from pathlib import Path

from biopipen.ns.delim import (
    RowsBinder as RowsBinder_,
    SampleInfo as SampleInfo_,
)
from biopipen.core.testing import get_pipeline


class RowsBinder(RowsBinder_):
    ...


class RowsBinderWithFilenames(RowsBinder_):
    envs = {"filenames": ["A", "B"]}


class RowsBinderWithFilenameStrings(RowsBinder_):
    envs = {"filenames": "A, B"}


class RowsBinderWithFilenameFuns(RowsBinder_):
    envs = {"filenames": "function(x) tools::file_path_sans_ext(basename(x))"}


class SampleInfo(SampleInfo_):
    requires = RowsBinderWithFilenameFuns
    envs = {
        "mutaters": {"Part": "paste0('Part_', Filename)"},
        "stats": {
            "Samples_Source": {
                "group": "Source",
            },
            "Samples_Subject": {
                "group": "Subject",
                "plot": "bar",
            },
            "Samples_Source_each_Subject": {
                "group": "Source",
                "each": "Subject",
            },
            "Samples_Subject_each_Source": {
                "group": "Subject",
                "plot": "bar",
                "each": "Source",
            },
            "Subjects_per_Filename": {
                "group": "Filename",
                "plot": "bar",
                "subset": "!duplicated(Subject)",
            },
            "Score": {
                "on": "Score",
                "subset": "!duplicated(Subject)"
            },
            "Score_Source": {
                "on": "Score",
                "group": "Source",
            },
            "Score_Source_each_Subject": {
                "on": "Score",
                "group": "Source",
                "each": "Subject",
            },
            "Score_violin": {
                "on": "Score",
                "plot": "violin",
                "subset": "!duplicated(Subject)"
            },
            "Score_Source_violin": {
                "on": "Score",
                "group": "Source",
                "plot": "violin",
            },
            "Score_Source_each_Subject_violin": {
                "on": "Score",
                "group": "Source",
                "each": "Subject",
                "plot": "violin",
            },
            "Score_vlnbox": {
                "on": "Score",
                "plot": "violin+box",
            },
            "Score_Source_vlnbox": {
                "on": "Score",
                "group": "Source",
                "plot": "violin+box",
            },
            "Score_Source_each_Subject_vlnbox": {
                "on": "Score",
                "group": "Source",
                "each": "Subject",
                "plot": "violin+box",
            },
            "Score_hist": {
                "on": "Score",
                "plot": "hist",
            },
            "Score_Source_hist": {
                "on": "Score",
                "group": "Source",
                "plot": "hist",
            },
            "Score_Source_each_Subject_hist": {
                "on": "Score",
                "group": "Source",
                "each": "Subject",
                "plot": "hist",
            },
        }
    }


class SampleInfo2(SampleInfo_):
    requires = RowsBinderWithFilenameFuns
    envs = {"exclude_cols": "Score,Filename"}


def pipeline():
    return (
        # get_pipeline(__file__, enable_report=True)
        get_pipeline(__file__)
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
    # assert pipen._succeeded
    outfile = (
        pipen.procs[-1].workdir.joinpath("0", "output", "A_rbound.txt")
    )
    assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
