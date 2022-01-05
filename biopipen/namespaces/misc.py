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


class Glob2Dir(Proc):
    """Create symbolic links in output directory for the files given
    by the glob pattern"""
    input = "pattern:var"
    output = "outdir:dir:from_glob"
    script = """
        for infile in {{in.pattern}}; do
            ln -s $(realpath $infile) "{{out.outdir}}/$(basename $infile)";
        done
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


class Str2File(Proc):
    """Write the given string to a file

    Input:
        str: The string to write to file
        name: The name of the file
            If not given, use `envs.name`

    Output:
        outfile: The output file

    Envs:
        name: The name of the output file
    """
    input = "str, name"
    output = "outfile:file:{{in.name}}"
    lang = config.lang.python
    envs = {"name": None}
    script = "file://../scripts/misc/Str2File.py"
