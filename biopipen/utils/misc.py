from __future__ import annotations
from pathlib import Path

import os
import sys
import logging
from subprocess import Popen
from typing import List, Callable, Any
from biopipen.core.filters import dict_to_cli_args  # noqa: F401

logger = logging.getLogger("biopipen_job")
logger.setLevel(logging.DEBUG)
_handler = logging.StreamHandler(sys.stdout)
# Use same log format as in R
# {sprintf("%-7s", level)} [{format(time, "%Y-%m-%d %H:%M:%S")}] {msg}
# so the logs can be populated by pipen-poplog
_handler.setFormatter(
    logging.Formatter(
        "%(levelname)-7s [%(asctime)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
)
logger.addHandler(_handler)


def require_package(package: str, version: str | None = None) -> None:
    """Require a Python package to be installed with optional version check.

    The version specifier should follow the format used by pip, e.g., '>=1.2.3'.
    Multiple version specifiers can be separated by commas, e.g., '>=1.2.3,<2.0.0'.

    Args:
        package (str): The name of the package to check.
        version (str | None): The version specifier string.
    """
    import importlib
    from importlib.metadata import version as get_version
    from packaging.specifiers import SpecifierSet

    try:
        importlib.import_module(package)
    except ImportError:
        raise ImportError(f"Package '{package}' is required but not installed.")

    if version:
        installed_version = get_version(package)
        specifier = SpecifierSet(version)
        if installed_version not in specifier:
            raise ImportError(
                f"Package '{package}' version '{installed_version}' does not satisfy "
                f"the requirement '{package}{version}'."
            )


def run_command(
    cmd: str | List[Any],
    fg: bool = False,
    wait: bool = True,
    print_command: bool = True,
    print_command_handler: Callable = print,
    **kwargs,
) -> Popen | str:
    """Run a command.

    Args:
        cmd: A string or list of strings representing the command to run.
        fg: Whether to run the command in the foreground.
            Redirects stdout and stderr to the current process.
        wait: Whether to wait for the command to finish.
            The command will be waited for if `fg` is `True`.
        print_command: Whether to print the command before running it.
        print_command_handler: The function to use to print the command.
        kwargs: Keyword arguments to pass to `subprocess.Popen`.

    Returns:
        The `Popen` object, or str when `stdout` is `RETURN` or `return`.
    """
    import shlex
    from subprocess import PIPE, STDOUT

    if isinstance(cmd, list):
        cmd = [str(c) for c in cmd]

    if print_command:
        print_command_handler("RUNNING COMMAND:")
        if isinstance(cmd, str):
            print_command_handler(f"  {cmd}\n")
        else:
            print_command_handler(f"  {shlex.join(cmd)}\n")
        # flush the output if print_command_handler is print
        if print_command_handler is print:
            sys.stdout.flush()

    if isinstance(cmd, str):
        kwargs["shell"] = True

    if kwargs.get("stdin") is True:
        kwargs["stdin"] = PIPE

    return_stdout = False
    stdout_file = None
    if kwargs.get("stdout") is True:
        kwargs["stdout"] = PIPE
    elif kwargs.get("stdout") in ("RETURN", "return"):
        kwargs["stdout"] = PIPE
        return_stdout = True
    elif isinstance(kwargs.get("stdout"), (str, Path)):
        if isinstance(kwargs["stdout"], str):
            kwargs["stdout"] = Path(kwargs["stdout"])
        stdout_file = kwargs["stdout"].open("w")
        kwargs["stdout"] = stdout_file
        kwargs["close_fds"] = True

    if kwargs.get("stderr") is True:
        kwargs["stderr"] = PIPE
    elif kwargs.get("stderr") in ("STDOUT", "stdout"):
        kwargs["stderr"] = STDOUT

    if fg:
        if kwargs.get("stdout") or kwargs.get("stderr"):
            raise ValueError(
                "Cannot redirect stdout or stderr when running in foreground"
            )
        kwargs["stdout"] = sys.stdout
        kwargs["stderr"] = sys.stderr
        kwargs["universal_newlines"] = True

    if "env" in kwargs:
        kwargs["env"] = {**os.environ, **kwargs["env"]}

    try:
        p = Popen(cmd, **kwargs)
    except Exception as e:
        raise RuntimeError(
            f"Failed to run command: {e}\n"
            f"Command (list): {cmd}\n"
            f"Command (str): {shlex.join(cmd)}"
        )

    if fg or wait or return_stdout:
        rc = p.wait()
        if rc != 0:
            if stdout_file:
                stdout_file.close()
            if return_stdout and p.stdout:
                p.stdout.close()
            raise RuntimeError(
                f"Failed to run command: rc={rc}\n"
                f"Command (list): {cmd}\n"
                f"Command (str): {shlex.join(cmd)}"
            )

        if return_stdout:
            try:
                return p.stdout.read().decode()  # type: ignore
            finally:
                p.stdout.close()  # type: ignore

        if stdout_file:
            stdout_file.close()

        return p

    return p
