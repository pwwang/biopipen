"""Script for bedtools.pBedGetfasta"""
# pylint: disable=undefined-variable,unused-import,invalid-name,not-a-mapping
# pylint: disable=unsubscriptable-object

from diot import Diot
from bioprocs.utils import shell2 as shell

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
bedtools = {{args.bedtools | quote}}
params = {{args.params | repr}}
ref = {{args.ref | quote}}

shell.load_config(bedtools=bedtools)

params.fi = ref
params.bed = infile
# params._out = outfile
# params._debug = True
if params['name+']:
    params.name = False

shell.bedtools.getfasta(**params).r > outfile
