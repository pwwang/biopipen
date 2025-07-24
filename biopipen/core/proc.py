"""Provides a base class for the processes to subclass"""
from __future__ import annotations

from diot import Diot  # type: ignore
from liquid.defaults import SEARCH_PATHS
from pipen import Proc as PipenProc  # type: ignore
from pipen_filters.filters import FILTERS

from .filters import filtermanager
from .defaults import BIOPIPEN_DIR, REPORT_DIR


def _repr(x):
    if isinstance(x, Diot):
        return repr(x.to_dict())
    return repr(x)


filtermanager.register("repr")(_repr)


class Proc(PipenProc):
    """Base class for all processes in biopipen to subclass"""

    template_opts = {
        "globals": {**FILTERS, "biopipen_dir": str(BIOPIPEN_DIR)},
        "filters": {**FILTERS, **filtermanager.filters},
        "search_paths": SEARCH_PATHS + [str(REPORT_DIR)],  # type: ignore
    }

    plugin_opts = {
        "poplog_pattern": (
            r"^(?P<level>INFO|WARN|WARNING|CRITICAL|ERROR|DEBUG?)\s*"
            r"\[\d+-\d+-\d+ \d+:\d+:\d+\] (?P<message>.*)$"
        )
    }
