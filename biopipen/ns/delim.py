"""Tools to deal with csv/tsv files"""
from ..core.config import config
from ..core.proc import Proc


class RowsBinder(Proc):
    """Bind rows of input files

    Input:
        infiles: The input files to bind.
            The input files should have the same number of columns,
            and same delimiter.

    Output:
        outfile: The output file with rows bound

    Envs:
        sep: The separator of the input files
        header (flag): Whether the input files have header
        filenames: Whether to add filename as the last column.
            Either a string of an R function that starts with `function` or
            a list of names (or string separated by comma) to
            add for each input file.
            The R function takes the path of the input file as the only
            argument and should return a string. The string will be added as
            the last column of the output file.
        filenames_col: The column name for the `filenames` columns
    """
    input = "infiles:files"
    output = (
        "outfile:file:"
        "{{in.infiles | first | stem}}_rbound{{in.infiles | first | ext}}"
    )
    envs = {
        "sep": "\t",
        "header": True,
        "filenames": None,
        "filenames_col": "Filename",
    }
    lang = config.lang.rscript
    script = "file://../scripts/delim/RowsBinder.R"
