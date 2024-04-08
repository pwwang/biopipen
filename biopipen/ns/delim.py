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


class SampleInfo(Proc):
    """List sample information and perform statistics

    Input:
        infile: The input file to list sample information
            The input file should be a csv/tsv file with header

    Output:
        outfile: The output file with sample information, with mutated columns
            if `envs.save_mutated` is True.

    Envs:
        sep: The separator of the input file.
        mutaters (type=json): A dict of mutaters to mutate the data frame.
            The key is the column name and the value is the R expression
            to mutate the column. The dict will be transformed to a list in R
            and passed to `dplyr::mutate`.
            You may also use `paired()` to identify paired samples. The function
            takes following arguments:
            * `df`: The data frame. Use `.` if the function is called in
                a dplyr pipe.
            * `id_col`: The column name in `df` for the ids to be returned in
                the final output.
            * `compare_col`: The column name in `df` to compare the values for
                each id in `id_col`.
            * `idents`: The values in `compare_col` to compare. It could be
                either an an integer or a vector. If it is an integer,
                the number of values in `compare_col` must be the same as
                the integer for the `id` to be regarded as paired. If it is
                a vector, the values in `compare_col` must be the same
                as the values in `idents` for the `id` to be regarded as paired.
            * `uniq`: Whether to return unique ids or not. Default is `TRUE`.
                If `FALSE`, you can mutate the meta data frame with the
                returned ids. Non-paired ids will be `NA`.
        save_mutated (flag): Whether to save the mutated columns.
        exclude_cols: The columns to exclude in the table in the report.
            Could be a list or a string separated by comma.
        defaults (ns): The default parameters for `envs.stats`.
            - on: The column name in the data for the stats.
                Default is `Sample`. The column could be either continuous or not.
            - subset: An R expression to subset the data.
                If you want to keep the distinct records, you can use
                `!duplicated(<col>)`.
            - group: The column name in the data for the group ids.
                If not provided, all records will be regarded as one group.
            - na_group (flag): Whether to include `NA`s in the group.
            - each: The column in the data to split the analysis in different
                plots.
            - ncol (type=int): The number of columns in the plot when `each`
                is not `NULL`. Default is 2.
            - na_each (flag): Whether to include `NA`s in the `each` column.
            - plot: Type of plot. If `on` is continuous, it could be
                `boxplot` (default), `violin`, `violin+boxplot` or `histogram`.
                If `on` is not continuous, it could be `barplot` or
                `pie` (default).
            - devpars (ns): The device parameters for the plot.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
        stats (type=json): The statistics to perform.
            The keys are the case names and the values are the parameters
            inheirted from `envs.defaults`.
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    envs = {
        "sep": "\t",
        "mutaters": {},
        "save_mutated": False,
        "exclude_cols": None,
        "defaults": {
            "on": "Sample",
            # "distinct": None,
            "group": None,
            "na_group": False,
            "each": None,
            "ncol": 2,
            "na_each": False,
            "plot": None,
            "devpars": {"width": 800, "height": 600, "res": 100},
        },
        "stats": {},
    }
    lang = config.lang.rscript
    script = "file://../scripts/delim/SampleInfo.R"
    plugin_opts = {"report": "file://../reports/delim/SampleInfo.svelte"}
