"""Misc processes"""
from ..core.proc import Proc

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
