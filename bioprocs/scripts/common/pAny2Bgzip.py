from bioprocs.utils import gztype
from bioprocs.utils import shell2 as shell

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

gzt = gztype(infile)
if gzt == 'gzip':
    shell.gunzip(c=infile, _pipe=True) | shell.bgzip(c=True, _out=outfile)
elif gzt == 'bgzip':
    shell.ln_s(infile, outfile)
else:
    shell.bgzip(c=infile, _out=outfile)
