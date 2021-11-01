"""Provides the envs from configuration files"""

import sys
from typing import Any
from pathlib import Path
from tempfile import gettempdir

from simpleconf import Config as SimConfig

from .defaults import BIOPIPEN_DIR

DEFAULT_CONFIG_FILE = BIOPIPEN_DIR / "core" / "config.toml"
USER_CONFIG_FILE = Path("~").expanduser() / ".biopipen.toml"
PROJ_CONFIG_FILE = Path(".") / ".biopipen.toml"


class Config(SimConfig):
    """Subclass the Config class to enables default value instead of
    KeyError for non-exist keys"""

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


config = Config(with_profile=False)

config_profiles = [
    {"path": {"tmpdir": gettempdir()}},
    DEFAULT_CONFIG_FILE,
    USER_CONFIG_FILE,
    PROJ_CONFIG_FILE,
]
# scan sys.argv to see if --config <config file> passed in
if "+config" in sys.argv:
    cindex = sys.argv.index("+config")
    config_profiles.append(sys.argv[cindex + 1])

config._load(*config_profiles)
