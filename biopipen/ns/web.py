"""Get data from the web"""

from ..core.proc import Proc
from ..core.config import config

class Download(Proc):
    """Download data from URLs

    Input:
        url: The URL to download data from

    Output:
        outfile: The file downloaded

    Envs:
        tool: Which tool to use to download the data
            wget, aria2c or python's urllib
        wget: Path to wget
        aria2c: Path to aria2c
        args: The arguments to pass to the tool
        ncores: The number of cores to use

    Requires:
        - name: wget
          message: Only required when envs.tool == "wget"
          check: |
            {{proc.envs.wget}} --version
        - name: aria2c
          message: Only required when envs.tool == "aria2c"
          check: |
            {{proc.envs.aria2c}} --version
    """
    input = "url"
    output = (
        "outfile:file:"
        "{{in.url | basename | replace: '%2E', '.' | slugify: separator='.'}}"
    )
    lang = config.lang.python
    envs = {
        "tool": "wget", # or aria2c, python
        "wget": config.exe.wget,
        "aria2c": config.exe.aria2c,
        "args": {},
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/web/Download.py"

class DownloadList(Proc):
    """Download data from URLs in a file

    Input:
        urlfile: The file containing the URLs to download data from

    Output:
        outdir: The directory containing the downloaded files

    Envs:
        tool: Which tool to use to download the data
            wget, aria2c or python's urllib
        wget: Path to wget
        aria2c: Path to aria2c
        args: The arguments to pass to the tool
        ncores: The number of cores to use

    Requires:
        - name: wget
          message: Only required when envs.tool == "wget"
          check: |
            {{proc.envs.wget}} --version
        - name: aria2c
          message: Only required when envs.tool == "aria2c"
          check: |
            {{proc.envs.aria2c}} --version
    """
    input = "urlfile:file"
    output = "outdir:dir:{{in.urlfile | stem}}.downloaded"
    lang = config.lang.python
    envs = {
        "tool": "wget", # or aria2c
        "wget": config.exe.wget,
        "aria2c": config.exe.aria2c,
        "args": {},
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/web/DownloadList.py"
