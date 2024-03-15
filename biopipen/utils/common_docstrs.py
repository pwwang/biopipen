"""Common docstrings for biopipen procs."""
import textwrap
from typing import Callable


def indent_docstr(docstr: str, indent: str) -> str:
    """Indent the docstring.

    Args:
        docstr: The docstring.
        indent: The indent.

    Returns:
        The indented docstring.
    """
    return textwrap.indent(docstr, indent).strip()


def format_placeholder(**kwargs) -> Callable[[type], type]:
    """A decorator to format a docstring placeholder.

    Args:
        **kwargs: The docstring placeholder.

    Returns:
        The decorated function.
    """

    def decorator(klass: type) -> type:
        klass.__doc__ = klass.__doc__ % kwargs
        return klass

    return decorator


MUTATE_HELPERS_CLONESIZE = """
There are also also 4 helper functions, `expanded`, `collapsed`, `emerged` and `vanished`,
which can be used to identify the expanded/collpased/emerged/vanished groups (i.e. TCR clones).
See also <https://pwwang.github.io/immunopipe/configurations/#mutater-helpers>.
For example, you can use
`{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}`
to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.
The values in this columns for other clones will be `NA`.
Those functions take following arguments:
* `df`: The metadata data frame. You can use the `.` to refer to it.
* `group.by`: The column name in metadata to group the cells.
* `idents`: The first group or both groups of cells to compare (value in `group.by` column). If only the first group is given, the rest of the cells (with non-NA in `group.by` column) will be used as the second group.
* `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).
* `each`: A column name (without quotes) in metadata to split the cells.
    Each comparison will be done for each value in this column (typically each patient or subject).
* `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).
* `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.
    If numeric column is given, the values should be the same for all cells in the same group.
    This will not be checked (only the first value is used).
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.
* `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.
* `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.
* `order`: The expression passed to `dplyr::arrange()` to order intermediate dataframe and get the ids in order accordingly.
    The intermediate dataframe includes the following columns:
    * `<id>`: The ids of clones (i.e. `CDR3.aa`).
    * `<each>`: The values in `each` column.
    * `ident_1`: The size of clones in the first group.
    * `ident_2`: The size of clones in the second group.
    * `.diff`: The difference between the sizes of clones in the first and second groups.
    * `.sum`: The sum of the sizes of clones in the first and second groups.
    * `.predicate`: Showing whether the clone is expanded/collapsed/emerged/vanished.
* `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.
* `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.

You can also use `top()` to get the top clones (i.e. the clones with the largest size) in each group.
For example, you can use
`{"Patient1_Top10_Clones": "top(subset = Patent == 'Patient1', uniq = FALSE)"}`
to create a new column in metadata named `Patient1_Top10_Clones`.
The values in this columns for other clones will be `NA`.
This function takes following arguments:
* `df`: The metadata data frame. You can use the `.` to refer to it.
* `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).
* `n`: The number of top clones to return. Default is `10`.
    If n < 1, it will be treated as the percentage of the size of the group.
    Specify `0` to get all clones.
* `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.
    If numeric column is given, the values should be the same for all cells in the same group.
    This will not be checked (only the first value is used).
    It is helpful to use `Clones` to use the raw clone size from TCR data, in case the cells are not completely mapped to RNA data.
    Also if you have `subset` set or `NA`s in `group.by` column, you should use `.n` to compare the number of cells in each group.
* `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).
* `each`: A column name (without quotes) in metadata to split the cells.
    Each comparison will be done for each value in this column (typically each patient or subject).
* `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.
* `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.
* `with_ties`: Whether to include ties (i.e. clones with the same size as the last clone) or not. Default is `FALSE`.
"""

ENVS_SECTION_EACH = """
The `section` is used to collect cases and put the results under the same directory and the same section in report.
When `each` for a case is specified, the `section` will be ignored and case name will be used as `section`.
The cases will be the expanded values in `each` column. When `prefix_each` is True, the column name specified by `each` will be prefixed to each value as directory name and expanded case name.
"""
