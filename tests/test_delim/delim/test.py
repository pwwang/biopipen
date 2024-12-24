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
                "plot_type": "pie",
                "x": "Source",
                "xlab": "",
                "ylab": "Number of Samples",
            },
            "Samples_Subject": {
                "plot_type": "pie",
                "x": "Subject",
                "xlab": "",
                "ylab": "Number of Samples",
            },
            "Samples_Source_each_Subject": {
                "plot_type": "bar",
                "x": "Source",
                "xlab": "",
                "ylab": "Number of Samples",
                "facet_by": "Subject",
            },
            "Samples_Subject_each_Source": {
                "plot_type": "bar",
                "x": "Subject",
                "xlab": "",
                "ylab": "Number of Samples",
                "facet_by": "Source",
            },
            "Subjects_per_Filename": {
                "plot_type": "bar",
                "x": "Filename",
                "subset": "!duplicated(Subject)",
            },
            "Score": {
                "plot_type": "density",
                "x": "Score",
                "subset": "!duplicated(Subject)",
            },
            "Score_Source": {
                "plot_type": "density",
                "x": "Score",
                "group_by": "Source",
            },
            "Score_Source_each_Subject": {
                "plot_type": "density",
                "x": "Score",
                "group_by": "Source",
                "facet_by": "Subject",
            },
            "Score_violin": {
                "plot_type": "violin",
                "y": "Score",
                "x": "Source",
                "subset": "!duplicated(Subject)",
            },
            "Score_Source_each_Subject_violin": {
                "plot_type": "violin",
                "y": "Score",
                "x": "Source",
                "facet_by": "Subject",
            },
            "Score_violin_box": {
                "plot_type": "violin",
                "y": "Score",
                "x": "Source",
                "add_box": True,
                "devpars": {"width": 800, "height": 600},
            },
            "Score_Source_each_Subject_violin_box": {
                "plot_type": "violin",
                "y": "Score",
                "x": "Source",
                "split_by": "Subject",
                "add_box": True,
            },
            "Score_hist": {
                "plot_type": "hist",
                "x": "Score",
            },
            "Score_Source_hist": {
                "plot_type": "hist",
                "x": "Score",
                "group_by": "Source",
            },
            "Score_Source_each_Subject_hist": {
                "plot_type": "hist",
                "x": "Score",
                "group_by": "Source",
                "facet_by": "Subject",
                "more_formats": "pdf",
                "save_code": True,
            },
        }
    }


class SampleInfo2(SampleInfo_):
    requires = RowsBinderWithFilenameFuns
    envs = {"exclude_cols": "Score,Filename"}


def pipeline():
    return (
        get_pipeline(__file__, enable_report=True)
        # get_pipeline(__file__)
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
