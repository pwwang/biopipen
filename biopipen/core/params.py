"""Provides the envs from configuration files"""

from typing import Any
from pathlib import Path

from simpleconf import Config

HERE = Path(__file__).parent.resolve()
DEFAULT_ENV_FILE = HERE / "params.toml"
USER_ENV_FILE = Path("~").expanduser() / ".biopipen.toml"
PROJ_ENV_FILE = Path(".") / ".biopipen.toml"

class Envs(Config):
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

params = Envs(with_profile=False)
params._load(DEFAULT_ENV_FILE, USER_ENV_FILE, PROJ_ENV_FILE)
