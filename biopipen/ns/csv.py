"""Tools to deal with csv/tsv files"""
from ..core.config import config
from ..core.proc import Proc


class BindRows(Proc):
    """Bind rows of input files"""
    input = "infiles:files"
    output = (
        "outfile:file:"
        "{{in.infiles | first | stem}}_rbind.{{in.infiles | first | ext}}"
    )
    envs = {
        "sep": "\t",
        "header": False,
        "helper_code": None,
        "add_filename": None,
    }
    lang = config.lang.python
    script = "file://../scripts/csv/BindRows.py"
