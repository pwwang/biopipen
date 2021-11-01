"""Misc processes"""
from ..core.proc import Proc

class File2Proc(Proc):
    """Accept a file and pass it down"""
    input = "infile:file"
    output = "outfile:file:{{in.infile | basename}}"
    script = "ln -s {{in.infile}} {{out.outfile}}"
