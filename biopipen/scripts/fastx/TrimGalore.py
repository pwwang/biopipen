"""Script for fastx.TrimGalore"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping
from pathlib import Path

from diot import Diot
from biopipen.utils import shell

fq1 = Path({{in.fq1 | repr}})
fq2 = Path({{in.fq2 | repr}})
fq1out = Path({{out.fq1 | repr}})
fq2out = Path({{out.fq2 | repr}})
trim_galore = {{args.trim_galore | repr}}
params = {{args.params | repr}}
outdir = Path({{job.outdir | repr}})

shell.load_config(trim_galore=trim_galore)

params.output_dir = outdir

shell.trim_galore(**params, _=[fq1, fq2]).fg()

# when paired = True, we only care about 1.R1_val_1.fq.gz not
# 1.R1_trimmed.fq.gz
if params.get('paired', False):
    for i, val in enumerate(outdir.glob('*_val_*.fq*')):
        if i == 0:
            if fq1out.is_file():
                fq1out.unlink()
            fq1out.symlink_to(val.resolve())
        else:
            if fq2out.is_file():
                fq2out.unlink()
            fq2out.symlink_to(val.resolve())
else:
    for i, trimmed in enumerate(outdir.glob('*_trimmed.fq*')):
        if i == 0:
            if fq1out.is_file():
                fq1out.symlink_to(trimmed)
        else:
            if fq2out.is_file():
                fq2out.symlink_to(trimmed)
