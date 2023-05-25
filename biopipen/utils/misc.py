from __future__ import annotations

import sys
from typing import Callable, List, Any


def exec_code(code, global_vars = None, local_vars = None, return_var = None):
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
    import shlex
    from subprocess import Popen, PIPE

    if isinstance(cmd, list):
        cmd = [str(c) for c in cmd]

    if print_command:
        print_command_handler("RUNNING COMMAND:")
        if isinstance(cmd, str):
            print_command_handler(f"  {cmd}")
        else:
            print_command_handler(f"  {shlex.join(cmd)}")

    if isinstance(cmd, str):
        kwargs["shell"] = True
    if kwargs.get("stdin") is True:
        kwargs["stdin"] = PIPE
    if kwargs.get("stdout") is True:
        kwargs["stdout"] = PIPE
    if kwargs.get("stderr") is True:
        kwargs["stderr"] = PIPE
    if fg:
        kwargs["stdout"] = sys.stdout
        kwargs["stderr"] = sys.stderr
        kwargs["universal_newlines"] = True
    p = Popen(cmd, **kwargs)
    if fg or wait:
        rc = p.wait()
        if rc != 0:
            raise RuntimeError(f"Failed to run command: {cmd}")
        return rc

    return p