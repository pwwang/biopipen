"""Get data from the web"""

from ..core.proc import Proc
from ..core.config import config

class Download(Proc):
    """Download data from URLs

    Input:
        url: The URL to download data from

    Output:
        outfile: The file downloaded
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
