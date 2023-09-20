"""Additional filters for pipen"""
from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any, List, Mapping

from argx import Namespace
from liquid.filters.manager import FilterManager

filtermanager = FilterManager()


@filtermanager.register
def dict_to_cli_args(
    dic: Mapping[str, Any],
    prefix: str | None = None,
    sep: str | None = " ",
    dup_key: bool = True,
    join: bool = False,
    start_key: str = "",
    end_key: str = "_",
    dashify: bool = False,
) -> str | List[str]:
    """Convert a python dict to a string of CLI arguments

    Args:
        dic: The dict to convert
        prefix: The prefix of the keys after conversion
            Defaults to `None`, mean `-` for short keys and `--` for long keys
        sep: The separator between key and value
            If `None`, using `" "` for short keys and `"="` for long keys
        dup_key: Whether to duplicate the key in cli arguments for list values
            When `True`, `{"a": [1, 2]}` will be converted to `"-a 1 -a 2"`
            When `False`, `{"a": [1, 2]}` will be converted to `"-a 1 2"`
            If `sep` is `None` or `=`, this must be True, otherwise an error
            will be raised
        join: Whether to join the arguments into a single string

    Returns:
        The converted string or list of strings
    """
    if sep in [None, "="] and not dup_key:
        raise ValueError("`dup_key` must be True when sep is `None` or `=`")

    starts = []
    ends = []
    out = []
    for k, v in dic.items():
        if k == start_key:
            container = starts
        elif k == end_key:
            container = ends
        else:
            container = out

        k = str(k)
        dashified_k = k.replace("_", "-") if dashify else k
        if v is None or v is False:
            continue

        if prefix is None:
            pref = "--" if len(k) > 1 else "-"
        else:
            pref = prefix

        if sep is None:
            s = "=" if len(k) > 1 else " "
        else:
            s = sep

        if v is True:
            # You can use {'-': True} to introduce a separator
            # like `--`
            if k in [start_key, end_key]:
                raise ValueError(
                    f"Cannot use `{start_key}` or `{end_key}` as key for True"
                )
            container.append(f"{pref}{dashified_k}")

        elif isinstance(v, (list, tuple)):
            for i, val in enumerate(v):
                if s == " ":
                    if (i == 0 or dup_key) and k not in [start_key, end_key]:
                        container.append(f"{pref}{dashified_k}")
                    container.append(str(val))
                else:
                    if (i == 0 or dup_key) and k not in [start_key, end_key]:
                        container.append(f"{pref}{dashified_k}{s}{val}")
                    else:
                        container.append(str(val))
        elif k in [start_key, end_key]:
            container.append(str(v))
        elif s == " ":
            container.append(f"{pref}{dashified_k}")
            container.append(str(v))
        else:
            container.append(f"{pref}{dashified_k}{s}{v}")

    out = starts + out + ends
    return shlex.join(out) if join else out


@filtermanager.register
def r(
    obj: Any,
    ignoreintkey: bool = True,
    todot: str = None,
    sortkeys: bool = False,
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
