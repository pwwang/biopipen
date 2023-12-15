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
For example, you can use
`{"Patient1_Tumor_Collapsed_Clones": "expanded(., Source, 'Tumor', subset = Patent == 'Patient1', uniq = FALSE)"}`
to create a new column in metadata named `Patient1_Tumor_Collapsed_Clones`
with the collapsed clones in the tumor sample (compared to the normal sample) of patient 1.
The values in this columns for other clones will be `NA`.
Those functions take following arguments:
* `df`: The metadata data frame. You can use the `.` to refer to it.
* `group-by`: The column name in metadata to group the cells.
* `idents`: The first group or both groups of cells to compare (value in `group-by` column). If only the first group is given, the rest of the cells (with non-NA in `group-by` column) will be used as the second group.
* `subset`: An expression to subset the cells, will be passed to `dplyr::filter()`. Default is `TRUE` (no filtering).
* `each`: A column name (without quotes) in metadata to split the cells.
    Each comparison will be done for each value in this column.
* `id`: The column name in metadata for the group ids (i.e. `CDR3.aa`).
* `compare`: Either a (numeric) column name (i.e. `Clones`) in metadata to compare between groups, or `.n` to compare the number of cells in each group.
    If numeric column is given, the values should be the same for all cells in the same group.
    This will not be checked (only the first value is used).
* `uniq`: Whether to return unique ids or not. Default is `TRUE`. If `FALSE`, you can mutate the meta data frame with the returned ids. For example, `df |> mutate(expanded = expanded(...))`.
* `debug`: Return the data frame with intermediate columns instead of the ids. Default is `FALSE`.
* `order`: The order of the returned ids. It could be `sum` or `diff`, which is the sum or diff of the `compare` between idents.
    Two kinds of modifiers can be added, including `desc` and `abs`.
    For example, `sum,desc` means the sum of `compare` between idents in descending order.
    Default is `diff,abs,desc`. It only works when `uniq` is `TRUE`. If `uniq` is `FALSE`, the returned
    ids will be in the same order as in `df`.
* `include_emerged`: Whether to include the emerged group for `expanded` (only works for `expanded`). Default is `FALSE`.
* `include_vanished`: Whether to include the vanished group for `collapsed` (only works for `collapsed`). Default is `FALSE`.
"""
