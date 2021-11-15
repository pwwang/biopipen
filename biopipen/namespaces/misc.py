"""Misc processes"""
from ..core.proc import Proc
from ..core.config import config

class File2Proc(Proc):
    """Accept a file and pass it down with a symbolic link

    Input:
        infile: The input file

    Output:
        outfile: The output symbolic link to the input file
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    script = """
        rm -f {{out.outfile | quote}}
        ln -s {{in.infile | quote}} {{out.outfile | quote}}
    """


class Config2File(Proc):
    """Write a configurationn in string to a configuration file

    Requires python package `toml`

    Input:
        config: A string representation of configuration
        name: The name for output file.
            Will be `config` if not given

    Output:
        outfile: The output file with the configuration

    Envs:
        infmt: The input format. `json` or `toml`.
        outfmt: The output format. `json` or `toml`.
    """
    input = "config:var, name:var"
    output = (
        "outfile:file:{{(in.name or 'config') | slugify}}.{{envs.outfmt}}"
    )
    envs = {"infmt": "toml", "outfmt": "toml"}
    lang = config.lang.python
    script = "file://../scripts/misc/Config2File.py"
