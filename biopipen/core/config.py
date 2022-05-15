"""Provides the envs from configuration files"""

import sys
from typing import Any
from pathlib import Path
from tempfile import gettempdir

from simpleconf import Config

from .defaults import BIOPIPEN_DIR

DEFAULT_CONFIG_FILE = BIOPIPEN_DIR / "core" / "config.toml"
USER_CONFIG_FILE = Path("~").expanduser() / ".biopipen.toml"
PROJ_CONFIG_FILE = Path(".") / ".biopipen.toml"

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

config = Config.load(*config_profiles, ignore_nonexist=True)
