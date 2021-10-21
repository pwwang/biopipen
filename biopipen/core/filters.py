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
