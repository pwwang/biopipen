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
            The basename of the output file will be the same as the input file.
            The file name of each plot will be slugified from the case name.
            Each plot has 3 formats: pdf, png and code.zip, which contains the
            data and R code to reproduce the plot.

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
        exclude_cols (auto): The columns to exclude in the table in the report.
            Could be a list or a string separated by comma.
        defaults (ns): The default parameters for `envs.stats`.
            - plot_type: The type of the plot.
                See the supported plot types here:
                <https://pwwang.github.io/plotthis/reference/index.html>
                The plot_type should be lower case and the plot function used in
                `plotthis` should be used. The mapping from plot_type to the
                plot function is like `bar -> BarPlot`, `box -> BoxPlot`, etc.
            - more_formats (list): The additional formats to save the plot.
                By default, the plot will be saved in png, which is also used to
                display in the report. You can add more formats to save the plot.
                For example, `more_formats = ["pdf", "svg"]`.
            - save_code (flag): Whether to save the R code to reproduce the plot.
                The data used to plot will also be saved.
            - subset: An expression to subset the data frame before plotting.
                The expression should be a string of R expression that will be passed
                to `dplyr::filter`. For example, `subset = "Sample == 'A'"`.
            - section: The section name in the report.
                In case you want to group the plots in the report.
            - devpars (ns): The device parameters for the plot.
                - width (type=int): The width of the plot.
                - height (type=int): The height of the plot.
                - res (type=int): The resolution of the plot.
            - descr: The description of the plot, shown in the report.
            - <more>: You can add more parameters to the defaults.
                These parameters will be expanded to the `envs.stats` for each case,
                and passed to individual plot functions.
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
            "plot_type": "bar",
            "more_formats": [],
            "save_code": False,
            "subset": None,
            "section": None,
            "descr": None,
            "devpars": {"width": None, "height": None, "res": 100},
        },
        "stats": {},
    }
    lang = config.lang.rscript
    script = "file://../scripts/delim/SampleInfo.R"
    plugin_opts = {"report": "file://../reports/common.svelte"}
