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
        # in case of deadlink
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
            if [[ -e $infile ]]; then
                ln -s $(realpath $infile) "{{out.outdir}}/$(basename $infile)";
            fi
        done
    """


class Config2File(Proc):
    """Write a configurationn in string to a configuration file

    Requires python package `rtoml`

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


class Shell(Proc):
    """Run a shell command

    Input:
        infile: The input file

    Output:
        outfile: The output file

    Envs:
        cmd: The shell command to run
            Use `$infile` and `$outfile` to refer to input and output files
        outdir: Whether the `out.outfile` should be a directory.
            If so a directory will be created before running the command.
    """
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    envs = {"cmd": "", "outdir": False}
    script = """
        infile={{in.infile | quote}}
        outfile={{out.outfile | quote}}
        is_outdir={{envs.outdir | int}}
        cmd={{envs.cmd | quote}}
        if [[ -z "$cmd" ]]; then
            echo "No command given." 1>&2
            exit 1
        fi
        if [[ $is_outdir -eq 1 ]]; then
            mkdir -p "$outfile"
        fi
        eval "$cmd"
    """
