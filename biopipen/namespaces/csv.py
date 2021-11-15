"""Tools to deal with csv/tsv files"""


from ..core.proc import Proc


class BindRows(Proc):
    """Bind rows of input files"""
    input = "infiles:files"
    output = "outfile:file:{{in.infiles[0] | stem}}.bound{{in.infiles[0] | ext}}"
    script = """
        outfile={{out.outfile | quote}}
        head -n 1 {{in.infiles[0] | quote}} > $outfile
        {% for infile in in.infiles %}
        tail -n +2 {{infile | quote}} >> $outfile
        {% endfor %}
    """
