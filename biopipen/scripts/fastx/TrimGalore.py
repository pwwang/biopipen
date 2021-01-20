"""Script for fastx.TrimGalore"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping

from diot import Diot
from biopipen.utils import shell

fq1 = {{in.fq1 | repr}}
fq2 = {{in.fq2 | repr}}
fq1out = {{out.fq1 | repr}}
fq2out = {{out.fq2 | repr}}
trim_galore = {{args.trim_galore | repr}}
params = {{args.params | repr}}
outdir = {{job.outdir | repr}}

shell.load_config(trim_galore=trim_galore)

params.output_dir = outdir

shell.trim_galore(**params, _=[fq1, fq2]).fg()
