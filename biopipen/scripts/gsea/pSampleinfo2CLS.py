
"""Script for gsea.pSampleinfo2CLS"""
# pylint: disable=undefined-variable,unused-import,invalid-name
from diot import OrderedDiot
from bioprocs.utils.sampleinfo import SampleInfo2 as SampleInfo

infile = {{i.sifile | quote}}
outfile = {{o.outfile | quote}}
samples = SampleInfo(infile).get_samples(return_all=True)
sam_to_group = {}
groups = []
for samrow in samples:
    if samrow['Group'] not in groups:
        groups.append(samrow['Group'])
    sam_to_group[samrow[0]] = samrow['Group']

with open(outfile, "w") as f:
    f.write("%s %s 1\n" % (len(samples), len(groups)))
    f.write("# %s\n" % (' '.join(groups)))
    f.write(' '.join(str(sam_to_group[samrow[0]])
                     for samrow in samples) + '\n')
