"""Script for vcf.pVcfPolyPhen"""

# pylint: disable=unused-import,invalid-name
from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell

# pylint: disable=undefined-variable
infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
polyphen2_annotator = {{args.polyphen2_annotator | quote}}
polyphen2_db = {{args.polyphen2_db | quote}}
# pylint: enable=undefined-variable

shell.load_config(polyphen2_annotator=polyphen2_annotator)

# Usage: vcf-annotate-polyphen polyphen.whess.sqlite input.vcf output.vcf

params = Diot()
params._ = [polyphen2_db, infile, outfile]

# If we got errors such like:
# SyntaxError: One of the FILTER lines is malformed: ##FILTER=<ID=ISNOTMUT,Description="Set if true: FMT/GT[0] = \"mis\" | FMT/GT[0] = \"0/0\"">
# Patch this to site-packages/vcf/parser.py (version: 0.6.8)
# 94c94
# <             Description="(?P<desc>[^"]*)"
# ---
# >             Description="(?P<desc>.*)"

# And path this to site-packages/vap/cli.py
# 52c52
# <                     'chr{}'.format(v.CHROM),
# ---
# >                     'chr{}'.format(v.CHROM.replace('chr', '')),

shell.polyphen2_annotator(**params).fg()
