from __future__ import annotations
from pathlib import Path

import os
import sys
import logging
from typing import List
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


def exec_code(code, global_vars=None, local_vars=None, return_var=None):
    global_vars = global_vars or {}
    local_vars = local_vars or {}
    exec(code, global_vars, local_vars)

    if return_var is not None:
        return local_vars[return_var]

    return None


def run_command(
    cmd: str | List[str],
    fg: bool = False,
    wait: bool = True,
    print_command: bool = True,
    print_command_handler: callable = print,
    **kwargs,
):
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
    from subprocess import Popen, PIPE, STDOUT

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
    if kwargs.get("stdout") is True:
        kwargs["stdout"] = PIPE
    elif kwargs.get("stdout") in ("RETURN", "return"):
        kwargs["stdout"] = PIPE
        return_stdout = True
    elif isinstance(kwargs.get("stdout"), (str, Path)):
        if isinstance(kwargs["stdout"], str):
            kwargs["stdout"] = Path(kwargs["stdout"])
        kwargs["stdout"] = kwargs["stdout"].open("w")
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
        raise RuntimeError(f"Failed to run command: {e}")

    if fg or wait or return_stdout:
        rc = p.wait()
        if rc != 0:
            raise RuntimeError(f"Failed to run command: {cmd}")

        if return_stdout:
            return p.stdout.read().decode()

        return p

    return p
