"""Script for vcf.pVcfSift4G"""

# pylint: disable=unused-import,invalid-name
from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell

# pylint: disable=undefined-variable
infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
sift4g_annotator = {{args.sift4g_annotator | quote}}
sift4g_db = {{args.sift4g_db | quote}}
params = {{args.params | repr}}
# pylint: enable=undefined-variable

shell.load_config(sift4g_annotator=sift4g_annotator)


# pylint: disable=line-too-long

# To run the SIFT 4G Annotator on Linux or Mac via command line, type the following command into the terminal:
# java -jar <Path to SIFT4G_Annotator> -c -i <Path to input vcf file> -d <Path to SIFT4G database directory> -r <Path to your results folder> -t
# Note:To run the Annotator via command line "-c" is essential (see the commandline parameters in the table below). If "-t" option is not used SIFT 4G extracts annotator single transcript per variant.

# Command line Options:

# Option	Description
# -c	To run on command line
# -i	Path to your input variants file in VCF format
# -d	Path to SIFT database directory
# -r	Path to your output results folder
# -t	To extract annotations for multiple transcripts (Optional)

# pylint: enable=line-too-long

# -c is not necessary, because sift4g_annotator has enabled it by default

params.i = infile
params.d = sift4g_db
params.r = Path(outfile).parent
params.t = True

# pylint: disable=not-a-mapping
shell.sift4g_annotator(**params).fg()
