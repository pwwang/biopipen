"""Additional filters for pipen"""
from pathlib import Path
from typing import Any, Mapping

import cmdy
from argx import Namespace
from liquid.filters.manager import FilterManager

filtermanager = FilterManager()


@filtermanager.register
def dict_to_cli_args(dic: Mapping[str, Any]) -> str:
    """Convert a python dict to a string of CLI arguments

    Examples:
        >>> {"a": 1, "ab": 2}
        >>> "-a 1 --ab 2"
    """
    return cmdy.echo(dic).stdout.strip()


@filtermanager.register
def r(
    obj: Any,
    ignoreintkey: bool = True,
    todot: str = None,
    sortkeys: bool = True,
) -> str:
    """Convert a python object into R repr

    Examples:
        >>> True -> "TRUE"
        >>> None -> "NULL"
        >>> [1, 2] -> c(1, 2)
        >>> {"a": 1, "b": 2} -> list(a = 1, b = 2)

    Args:
        ignoreintkey: When keys of a dict are integers, whether we should
            ignore them. For example, when `True`, `{1: 1, 2: 2}` will be
            translated into `"list(1, 2)"`, but `"list(`1` = 1, `2` = 2)"`
            when `False`
        todot: If not None, the string will be converted to a dot
            For example, `todot="-"` will convert `"a-b"` to `"a.b"`
            Only applies to the keys of obj when it is a dict
        sortkeys: Whether to sort the keys of a dict.
            True by default, in case the order of keys matters, for example,
            it could affect whether a job is cached.
            But sometimes, you want to keep orginal order, for example,
            arguments passed the `dplyr::mutate` function. Because the later
            arguments can refer to the earlier ones.

    Returns:
        Then converted string representation of the object
    """
    if obj is True:
        return "TRUE"
    if obj is False:
        return "FALSE"
    if obj is None:
        return "NULL"
    if isinstance(obj, str):
        if obj.upper() in ["+INF", "INF"]:
            return "Inf"
        if obj.upper() == "-INF":
            return "-Inf"
        if obj.upper() == "TRUE":
            return "TRUE"
        if obj.upper() == "FALSE":
            return "FALSE"
        if obj.upper() == "NA" or obj.upper() == "NULL":
            return obj.upper()
        if obj.startswith("r:") or obj.startswith("R:"):
            return str(obj)[2:]
        return repr(str(obj))
    if isinstance(obj, Path):
        return repr(str(obj))
    if isinstance(obj, (list, tuple, set)):
        if any(isinstance(i, dict) for i in obj):
            # c(list(a=1), list(b=2)) will be combined as list(a=1, b=2)
            # but we want list(list(a=1), list(b=2))
            wrapper = "list"
        else:
            wrapper = "c"
        return "{}({})".format(
            wrapper,
            ",".join([r(i, ignoreintkey, todot, sortkeys) for i in obj]),
        )
    if isinstance(obj, dict):
        # list allow repeated names
        return "list({})".format(
            ",".join(
                [
                    "`{0}`={1}".format(k, r(v, ignoreintkey, todot, sortkeys))
                    if isinstance(k, int) and not ignoreintkey
                    else r(v, ignoreintkey, todot, sortkeys)
                    if isinstance(k, int) and ignoreintkey
                    else "`{0}`={1}".format(
                        (
                            str(k)
                            if not todot
                            else str(k).replace(todot, ".")
                        ).split("#")[0],
                        r(v, ignoreintkey, todot, sortkeys),
                    )
                    for k, v in (
                        sorted(obj.items()) if sortkeys else obj.items()
                    )
                ]
            )
        )
    if isinstance(obj, Namespace):
        return r(vars(obj), ignoreintkey, todot, sortkeys)
    return repr(obj)
