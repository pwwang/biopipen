"""Script for pCnvkitFlatRef"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping,unsupported-assignment-operation
# pylint: disable=unsubscriptable-object

from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell

cnvkit   = {{args.cnvkit | quote}}
infile   = {{i.tgfile | quote}}
atgfile  = {{i.atgfile | quote}}
ref      = {{args.ref | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
nthread  = 1

shell.load_config(cnvkit = dict(
    _exe = cnvkit,
    _env = dict(
        OPENBLAS_NUM_THREADS = str(nthread),
        OMP_NUM_THREADS      = str(nthread),
        NUMEXPR_NUM_THREADS  = str(nthread),
        MKL_NUM_THREADS      = str(nthread)
    ),
    _cwd = path.dirname(outfile)
))

params.o = outfile
params.f = ref
params.t = infile
params.a = atgfile
shell.cnvkit.reference(**params).fg
