"""Provides the envs from configuration files"""

from typing import Any
from pathlib import Path
from tempfile import gettempdir

from diot import Diot
from simpleconf import Config

from .defaults import BIOPIPEN_DIR

DEFAULT_CONFIG_FILE = BIOPIPEN_DIR / "core" / "config.toml"
USER_CONFIG_FILE = Path("~").expanduser() / ".biopipen.toml"
PROJ_CONFIG_FILE = Path(".") / ".biopipen.toml"


class ConfigItems(Diot):
    """Provides the envs from configuration files and defaults the
    non-existing values to None."""

    def __getattr__(self, name: str) -> Any:
        try:
            return super().__getattr__(name)
        except (KeyError, AttributeError):
            return None

    def __getitem__(self, name: str) -> Any:
        try:
            return super().__getitem__(name)
        except (KeyError, AttributeError):
            return None


config_profiles = [
    {"path": {"tmpdir": gettempdir()}},
    DEFAULT_CONFIG_FILE,
    USER_CONFIG_FILE,
    PROJ_CONFIG_FILE,
]

config = ConfigItems(Config.load(*config_profiles, ignore_nonexist=True))
