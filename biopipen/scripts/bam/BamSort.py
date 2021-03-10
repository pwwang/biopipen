"""Script for bam.PatternCNV"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping
from pathlib import Path, PosixPath

from diot import Diot
from biopipen.utils import shell

bamfile = {{in.bamfile | repr}}
outfile = Path({{out.outfile | repr}})
tool = {{args.tool | repr}}
samtools = {{args.samtools | repr}}
ref = {{args.ref | repr}}
ncores = {{args.ncores | repr}}
by = {{args.by | repr}}
params = {{args.params | repr}}

shell.load_config(samtools=samtools)

def tool_samtools():
    params.n = by == 'name'

    shell.samtools.sort(
        _=bamfile,
        o=outfile,
        reference=ref,
        threads=ncores,
        **params
    ).fg()

globals()[f'tool_{tool}']()
