"""Provides a base class for the processes to subclass"""

from pathlib import Path

from pipen import Proc as PipenProc

from .filters import filtermanager

class Proc(PipenProc):
    """Base class for all processes in biopipen to subclass"""
    template_opts = {
        "globals": {
            "biopipen_dir": str(Path(__file__).parent.parent),
        },
        "filters": filtermanager.filters.copy(),
    }
