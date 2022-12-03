"""Provides a base class for the processes to subclass"""
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Mapping

from liquid.defaults import SEARCH_PATHS
from diot import Diot
from pipen import Pipen, Proc as PipenProc

from .config import config
from .filters import filtermanager
from .defaults import BIOPIPEN_DIR, REPORT_DIR


class Proc(PipenProc):
    """Base class for all processes in biopipen to subclass"""

    template_opts = {
        "globals": {
            "biopipen_dir": str(BIOPIPEN_DIR),
        },
        "filters": filtermanager.filters.copy(),
        "search_paths": SEARCH_PATHS + [str(REPORT_DIR)],
    }


class Pipeline(ABC):
    """The pipeline base class"""

    __slots__ = ("options", "starts", "ends", "procs")

    _PIPELING = None
    defaults = Diot()

    def __new__(cls, *args, **kwargs):
        if cls._PIPELING is None:
            cls._PIPELING = super().__new__(cls)
        return cls._PIPELING

    def __init__(self, options: Mapping[str, Any] | None = None):
        self.options = (
            Diot(self.__class__.defaults)
            | (config.get("pipeline", {}).get("cnvkit_pipeline", {}))
            | (options or {})
        )
        self.starts = []
        self.ends = []
        self.procs = Diot()

        self.build()

        if not self.starts:
            raise ValueError(
                "No start processes found. "
                "Did you forget to add them in build()?"
            )

    @abstractmethod
    def build(self) -> None:
        """Build the pipeline"""
        raise NotImplementedError

    def run(self, data=None, **kwargs) -> Pipen:
        """Run the pipeline

        Args:
            data: The data to run the pipeline with
            kwargs: The keyword arguments to pass to Pipen

        Returns:
            The Pipen instance
        """

        pipe = Pipen(**kwargs).set_start(self.starts)
        if data is not None:
            pipe.set_data(data)
        pipe.run()
        return pipe
