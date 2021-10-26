from pathlib import Path
from typing import Any, Mapping

import cmdy
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
def r(obj: Any, ignoreintkey: bool = True) -> str:
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

    Returns:
        Then converted string representation of the object
    """
    if obj is True:
        return 'TRUE'
    if obj is False:
        return 'FALSE'
    if obj is None:
        return 'NULL'
    if isinstance(obj, str):
        if obj.upper() in ['+INF', 'INF']:
            return 'Inf'
        if obj.upper() == '-INF':
            return '-Inf'
        if obj.upper() == 'TRUE':
            return 'TRUE'
        if obj.upper() == 'FALSE':
            return 'FALSE'
        if obj.upper() == 'NA' or obj.upper() == 'NULL':
            return obj.upper()
        if obj.startswith('r:') or obj.startswith('R:'):
            return str(obj)[2:]
        return repr(str(obj))
    if isinstance(obj, Path):
        return repr(str(obj))
    if isinstance(obj, (list, tuple, set)):
        return 'c({})'.format(','.join([r(i) for i in obj]))
    if isinstance(obj, dict):
        # list allow repeated names
        return 'list({})'.format(','.join([
            '`{0}`={1}'.format(
                k,
                r(v)) if isinstance(k, int) and not ignoreintkey else \
                r(v) if isinstance(k, int) and ignoreintkey else \
                '`{0}`={1}'.format(str(k).split('#')[0], r(v))
            for k, v in sorted(obj.items())]))
    return repr(obj)
