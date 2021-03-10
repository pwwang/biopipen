"""Script for cnvkit.CNVkitBatch"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping
from pathlib import Path, PosixPath

from diot import Diot
from biopipen.utils import shell

tumors = {{in.tumors | repr}}
normals = {{in.normals | repr}}
outdir = Path({{out.outdir | repr}})
cnvkit = {{args.cnvkit | repr}}
ncores = {{args.ncores | repr}}
ref = {{args.ref | repr}}
params = {{args.params | repr}}

shell.mkdir(p=outdir)

shell.load_config(cnvkit=cnvkit)

# generate access if not provided
# if params.get('method') != 'wgs':
#     params.access = params.access or 5000
#     if isinstance(params.access, int):
#         shell.cnvkit.access(ref, s=params.access, o=outdir / 'access.bed').fg()
#         params.access = outdir / 'access.bed'

params.normal = normals or True
params.fasta = ref
params['output-reference'] = outdir / 'reference.cnn'
params['output-dir'] = outdir
params['p'] = ncores

if params.targets is None:
    del params['targets']

if params.access is None:
    del params['access']

shell.cnvkit.batch(*tumors, **params).fg()
