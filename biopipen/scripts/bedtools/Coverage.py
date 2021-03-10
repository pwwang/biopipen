"""Script for bedtools.Coverage"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping

from diot import Diot
from biopipen.utils import shell

afile = {{in.afile | repr}}
bfile = {{in.bfile | repr}}
outfile = {{out.outfile | repr}}
bedtools = {{args.bedtools | repr}}
params = {{args.params | repr}}

shell.load_config(bedtools=bedtools)

params.a = afile
params.b = bfile
shell.bedtools.coverage(**params).r() > outfile
