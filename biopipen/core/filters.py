"""Additional filters for pipen"""
from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any, List, Mapping

from argx import Namespace
from liquid.filters.manager import FilterManager
from pipen_report.filters import register_component, render_ui, _tag

filtermanager = FilterManager()


@filtermanager.register
def dict_to_cli_args(
    dic: Mapping[str, Any],
    exclude: List[str] = None,
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
        exclude: The keys to exclude
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
        start_key: The key to start the arguments
            This is useful when you want to put some arguments at the beginning
            of the command line
        end_key: The key to end the arguments
            This is useful when you want to put some arguments at the end
            of the command line
        dashify: Whether to replace `_` with `-` in the keys

    Returns:
        The converted string or list of strings
    """
    if sep in [None, "="] and not dup_key:
        raise ValueError("`dup_key` must be True when sep is `None` or `=`")

    if exclude:
        dic = {k: v for k, v in dic.items() if k not in exclude}

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
    skip: int = 0,
    _i: int = 0,
) -> str:
    """Convert a python object into R repr

    Examples:
        >>> True -> "TRUE"
        >>> None -> "NULL"
        >>> [1, 2] -> c(1, 2)
        >>> {"a": 1, "b": 2} -> list(a = 1, b = 2)

    Args:
        obj: The object to convert
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
        skip: Levels to skip for `todot`. For example, `skip=1` will skip
            the first level of the keys. When `todot` is `"-"`, `skip=1` will
            convert `{"a-b": {"c-d": 1}}` to ``list(`a-b` = list(`c.d` = 1))``
        _i: Current level of the keys. Used internally

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
            ", ".join(
                [r(i, ignoreintkey, todot, sortkeys, skip, _i + 1) for i in obj]
            ),
        )
    if isinstance(obj, dict):
        # list allow repeated names
        items = []
        keys = obj.keys()
        if sortkeys:
            keys = sorted(keys)
        for k in keys:
            v = obj[k]
            if isinstance(k, int) and not ignoreintkey:
                item = (
                    f"`{k}`={r(v, ignoreintkey, todot, sortkeys, skip, _i + 1)}"
                )
            elif isinstance(k, int) and ignoreintkey:
                item = r(v, ignoreintkey, todot, sortkeys, skip, _i + 1)
            else:
                key = str(k)
                if todot and _i >= skip:
                    key = key.replace(todot, ".")
                item = (
                    f"`{key}`="
                    f"{r(v, ignoreintkey, todot, sortkeys, skip, _i + 1)}"
                )
            items.append(item)

        return f"list({', '.join(items)})"

    if isinstance(obj, Namespace):
        return r(vars(obj), ignoreintkey, todot, sortkeys, skip, _i)

    return repr(obj)


@register_component("fgsea")
def _render_fgsea(
    cont: Mapping[str, Any],
    job: Mapping[str, Any],
    level: int,
    na_arg: str = "10",
) -> str:
    """Render fgsea report"""
    # cont["dir"] is required
    n_pathways = int(na_arg)
    pathways = []
    with Path(cont["dir"]).joinpath("fgsea.txt").open() as f:
        next(f)  # skip header
        for line in f:
            items = line.strip().split("\t")
            pathways.append((items[0], items[-1]))
            if len(pathways) >= n_pathways:
                break

    components = [
        # Summary
        {
            "title": "Enrichment Analysis Summary",
            "ui": "tabs",
            "contents": [
                {
                    "title": "Plot",
                    "ui": "flat",
                    "contents": [
                        {
                            "kind": "image",
                            "src": str(Path(cont["dir"]).joinpath("gsea_table.png")),
                        }
                    ],
                },
                {
                    "title": "Table",
                    "ui": "flat",
                    "contents": [
                        {
                            "kind": "table",
                            "src": str(Path(cont["dir"]).joinpath("fgsea.txt")),
                            "data": {"excluded": {"slug"}},
                        }
                    ],
                },
            ]
        },
        # Pathways
        {
            "title": f"Enriched Pathways (Top {n_pathways})",
            "ui": "table_of_images",
            "contents": [
                {
                    "src": str(Path(cont["dir"]) / f"fgsea_{slug}.png"),
                    "title": pw,
                }
                for pw, slug in pathways
            ]
        },
    ]

    return render_ui(components, "accordion", job, level)


@register_component("pdf")
def _render_pdf(
    cont: Mapping[str, Any],
    job: Mapping[str, Any],
    level: int,
) -> str:
    """Render pdf report"""
    # cont["src"] is required
    height = cont.get("height", "600")
    return _tag(
        "embed",
        src=str(cont["src"]),
        type="application/pdf",
        width="100%",
        height=height,
    )


@register_component("gsea")
def _render_gsea(
    cont: Mapping[str, Any],
    job: Mapping[str, Any],
    level: int,
) -> str:
    """Render gsea report"""
    # cont["dir"] is required
    raise NotImplementedError()


@register_component("enrichr")
def _render_enrichr(
    cont: Mapping[str, Any],
    job: Mapping[str, Any],
    level: int,
) -> str:
    """Render enrichr report"""
    # cont["dir"] is required
    dbs = [sumfile.stem[8:] for sumfile in Path(cont["dir"]).glob("Enrichr-*.txt")]
    components = []

    for db in dbs:
        enrichr_plot = Path(cont["dir"]).joinpath(f"Enrichr-{db}.png")
        if enrichr_plot.exists():
            components.append(
                {
                    "title": db,
                    "ui": "tabs",
                    "contents": [
                        {
                            "title": "Plot",
                            "ui": "flat",
                            "contents": [
                                {
                                    "kind": "image",
                                    "src": str(
                                        Path(cont["dir"]).joinpath(f"Enrichr-{db}.png")
                                    ),
                                }
                            ],
                        },
                        {
                            "title": "Table",
                            "ui": "flat",
                            "contents": [
                                {
                                    "kind": "table",
                                    "src": str(
                                        Path(cont["dir"]).joinpath(f"Enrichr-{db}.txt")
                                    ),
                                }
                            ],
                        },
                    ],
                }
            )
        else:
            components.append(
                {
                    "title": db,
                    "ui": "tabs",
                    "contents": [
                        {
                            "title": "Error",
                            "ui": "flat",
                            "contents": [
                                {
                                    "kind": "error",
                                    "content": "No enriched terms found.",
                                }
                            ],
                        },
                    ],
                }
            )

    return render_ui(components, "accordion", job, level)
