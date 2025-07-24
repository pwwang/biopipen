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
    lang = config.lang.bash
    script = """
        # in case of deadlink
        rm -f {{out.outfile | quote}}
        if [[ ! -e {{in.infile | quote}} ]]; then
            echo "File {{in.infile | quote}} does not exist." 1>&2
            exit 1
        fi
        ln -s {{in.infile | quote}} {{out.outfile | quote}}
    """


class Glob2Dir(Proc):
    """Create symbolic links in output directory for the files given
    by the glob pattern"""
    input = "pattern:var"
    output = "outdir:dir:from_glob"
    lang = config.lang.bash
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
    output = "outfile:file:{{(in.name or 'config') | slugify}}.{{envs.outfmt}}"
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
    output = "outfile:file:{{in.name | default: 'unnamed.txt'}}"
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
    lang = config.lang.bash
    script = "file://../scripts/misc/Shell.sh"


class Plot(Proc):
    """Plot given data using plotthis package in R

    Input:
        datafile: The input data file in RDS or qs/qs2 format.
            If it is not in RDS nor qs/qs2 format, read.table will be used
            to read the data file with the options provided by `envs.read_opts`.

    Output:
        plotfile: The output plot file in PNG format

    envs:
        fn: The plot function to use. Required.
        devpars (ns): The device parameters for the plot.
            - width: The width of the plot in pixels.
            - height: The height of the plot in pixels.
            - res: The resolution of the plot in DPI.
        more_formats: The additional formats to save the plot in other than PNG.
            The file will be saved in the same directory as the plotfile.
        save_code: Whether to save the R code used for plotting.
        read_opts: Options to read the data file.
            If the data file is not in RDS nor qs/qs2 format, these options
            will be passed to `read.table`.
        <more>: Additional parameters to the plot function.
    """
    input = "datafile:file"
    output = "plotfile:file:{{in.datafile | stem}}.png"
    envs = {
        "fn": None,
        "devpars": {"res": 100},
        "more_formats": [],
        "save_code": False,
        "read_opts": {},
    }
    lang = config.lang.rscript
    script = "file://../scripts/misc/Plot.R"
